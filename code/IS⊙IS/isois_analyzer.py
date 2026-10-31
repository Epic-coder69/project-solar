#!/usr/bin/env python3
"""Parker Solar Probe IS⊙IS unified analyzer.

Contains actual survey-rate, fine-event and multiresolution coverage analysis
algorithms, not subprocess proxies. Nominal timestamp overlap is NOT proof of
particle quality, statistically significant counts or physical acceleration.
Previous output files and formats are preserved.
"""
from __future__ import annotations

from isois_reader import (
    _atomic_csv, _atomic_text, _sha256, _utc,
    build_event_windows, validate_fields_events,
    EPIHI_DEFAULT_PRODUCTS as DEFAULT_PRODUCTS,
)

# === SUMMARY-RATE SURVEY ANALYSIS (from particle_data_analysis.py) ===
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

READER_COLUMNS = {
    'timestamp_utc', 'integration_start_utc', 'integration_end_utc',
    'variable', 'value', 'units', 'valid_measurement', 'product', 'analysis_grade',
}
EVENT_COLUMNS = {'event_id', 'event_start_utc', 'event_end_utc'}


def _survey_utc(series):
    return pd.to_datetime(series, utc=True, errors='raise')


def validate_particles(particles):
    missing = READER_COLUMNS - set(particles.columns)
    if missing:
        raise ValueError(f'Missing reader columns: {sorted(missing)}')
    if particles.empty:
        raise ValueError('No particle samples')
    if set(particles['analysis_grade'].dropna()) != {'survey_only'}:
        raise ValueError('This exploratory analysis only accepts labelled survey_only input')
    if set(particles['product'].dropna()) != {'PSP_ISOIS_L2-SUMMARY'}:
        raise ValueError('Unrecognized particle product; cannot silently mix instruments')
    frame = particles.copy()
    for col in ('timestamp_utc', 'integration_start_utc', 'integration_end_utc'):
        frame[col] = _survey_utc(frame[col])
    if frame[['timestamp_utc', 'integration_start_utc', 'integration_end_utc']].isna().any().any():
        raise ValueError('Missing timestamps')
    valid = frame['valid_measurement'].astype(str).str.lower()
    if not valid.isin(['true', 'false']).all():
        raise ValueError('Invalid boolean valid_measurement values')
    frame['valid_measurement'] = valid.eq('true')
    frame['value'] = pd.to_numeric(frame['value'], errors='coerce')
    if ((frame.integration_start_utc >= frame.integration_end_utc) |
        (frame.timestamp_utc < frame.integration_start_utc) |
        (frame.timestamp_utc > frame.integration_end_utc)).any():
        raise ValueError('Invalid or non-midpoint sample integration intervals')
    if (frame.loc[frame.valid_measurement, 'value'].isna().any() or
        (~np.isfinite(frame.loc[frame.valid_measurement, 'value'])).any() or
        (frame.loc[frame.valid_measurement, 'value'] < 0).any()):
        raise ValueError('Valid samples have nonfinite or negative values')
    if (frame.loc[~frame.valid_measurement, 'value'].notna()).any():
        raise ValueError('Invalid samples must have missing values, not zero')
    if frame.duplicated(['variable', 'timestamp_utc']).any():
        raise ValueError('Duplicate particle epochs')
    for variable, group in frame.groupby('variable'):
        if len(group.units.dropna().unique()) != 1:
            raise ValueError(f'{variable}: inconsistent units between days')
        g = group.sort_values('integration_start_utc')
        if (g.integration_start_utc.iloc[1:].to_numpy() <
            g.integration_end_utc.iloc[:-1].to_numpy()).any():
            raise ValueError(f'{variable}: overlapping integrations (would double-count exposure)')
    return frame


def summarize_exposure(group, left, right, *, min_coverage=0.7, min_samples=3):
    """Weight by actual overlap and cover only observed valid seconds.

    The coverage uses the union of valid integration intervals; overlap is
    forbidden by validation, so summing durations is valid here.
    """
    if not left < right:
        raise ValueError('Nonpositive target interval')
    first_idx = int(group.integration_end_utc.searchsorted(left, side='right'))
    last_idx = int(group.integration_start_utc.searchsorted(right, side='left'))
    group = group.iloc[first_idx:last_idx]
    a = group.integration_start_utc.clip(lower=left)
    b = group.integration_end_utc.clip(upper=right)
    overlap = (b - a).dt.total_seconds().clip(lower=0).to_numpy(dtype=float)
    values = group.value.to_numpy(dtype=float)
    valid = group.valid_measurement.to_numpy(dtype=bool) & (overlap > 0) & np.isfinite(values)
    seconds = float(overlap[valid].sum())
    duration = (right - left).total_seconds()
    coverage = min(seconds / duration, 1.0)
    n = int(valid.sum())
    ok = n >= min_samples and coverage >= min_coverage and seconds > 0
    weighted = float(np.average(values[valid], weights=overlap[valid])) if ok else np.nan
    return {'mean_rate': weighted, 'coverage_fraction': coverage,
            'valid_samples': n, 'valid_seconds': seconds, 'usable': ok}


def characterize(events, particles, *, context_minutes=30, min_coverage=0.7,
                 min_samples=3):
    if not EVENT_COLUMNS.issubset(events.columns):
        raise ValueError(f'FIELDS catalog missing {sorted(EVENT_COLUMNS - set(events.columns))}')
    if events.empty or events.event_id.isna().any() or events.event_id.duplicated().any():
        raise ValueError('Event catalog is empty or contains missing/duplicate IDs')
    if not 0 <= min_coverage <= 1 or min_samples < 1 or context_minutes <= 0:
        raise ValueError('Invalid analysis thresholds')
    obs = validate_particles(particles)
    out = []
    for variable, group in obs.groupby('variable'):
        group = group.sort_values('integration_start_utc').reset_index(drop=True)
        units = str(group.units.iloc[0])
        for event in events.to_dict('records'):
            first = pd.Timestamp(event['event_start_utc'])
            last = pd.Timestamp(event['event_end_utc'])
            first = first.tz_localize('UTC') if first.tzinfo is None else first.tz_convert('UTC')
            last = last.tz_localize('UTC') if last.tzinfo is None else last.tz_convert('UTC')
            if last <= first:
                raise ValueError(f'Nonpositive FIELDS event interval: {event["event_id"]}')
            pad = pd.Timedelta(minutes=context_minutes)
            windows = {
                'pre': (first - pad, first),
                'during': (first, last),
                'post': (last, last + pad),
            }
            record = {
                'event_id': event['event_id'],
                'event_start_utc': first.isoformat(),
                'event_end_utc': last.isoformat(),
                'variable': variable, 'units': units,
                'product': 'PSP_ISOIS_L2-SUMMARY', 'analysis_grade': 'survey_only',
                'context_minutes': context_minutes,
            }
            usable = {}
            for label, (left, right) in windows.items():
                stat = summarize_exposure(group, left, right,
                                          min_coverage=min_coverage, min_samples=min_samples)
                record[f'{label}_mean_rate'] = stat['mean_rate']
                record[f'{label}_coverage_fraction'] = stat['coverage_fraction']
                record[f'{label}_valid_samples'] = stat['valid_samples']
                record[f'{label}_valid_seconds'] = stat['valid_seconds']
                usable[label] = stat['usable']
            if not all(usable.values()):
                record['status'] = 'insufficient_particle_coverage'
                record['during_vs_pre_pct'] = np.nan
            elif record['pre_mean_rate'] <= 0:
                record['status'] = 'zero_pre_baseline'
                record['during_vs_pre_pct'] = np.nan
            else:
                record['status'] = 'survey_characterized'
                record['during_vs_pre_pct'] = (
                    100 * (record['during_mean_rate'] / record['pre_mean_rate'] - 1))
            out.append(record)
    return pd.DataFrame(out)


def _survey_self_test():
    start = pd.Timestamp('2018-11-06T12:00:00Z')
    rows = []
    for i in range(-40, 41):
        midpoint = start + pd.Timedelta(minutes=i)
        rows.append({
            'timestamp_utc': midpoint.isoformat(),
            'integration_start_utc': (midpoint - pd.Timedelta(seconds=30)).isoformat(),
            'integration_end_utc': (midpoint + pd.Timedelta(seconds=30)).isoformat(),
            'variable': 'A_H_Rate_TS', 'value': 20.0 if 0 <= i < 10 else 10.0,
            'valid_measurement': True, 'units': 'counts s!E-1!N',
            'product': 'PSP_ISOIS_L2-SUMMARY', 'analysis_grade': 'survey_only',
        })
    events = pd.DataFrame([{'event_id': 'PSP-MAG-test',
                            'event_start_utc': start.isoformat(),
                            'event_end_utc': (start + pd.Timedelta(minutes=10)).isoformat()}])
    p = pd.DataFrame(rows)
    row = characterize(events, p, min_samples=3).iloc[0]
    assert np.isclose(row.pre_mean_rate, (29.5 * 10.0 + .5 * 20.0) / 30.0)
    assert 18.0 < row.during_mean_rate < 20.0
    assert row.status == 'survey_characterized'
    assert 0.99 < row.pre_coverage_fraction <= 1.0
    p.loc[p['timestamp_utc'] == start.isoformat(), 'value'] = np.nan
    p.loc[p['timestamp_utc'] == start.isoformat(), 'valid_measurement'] = False
    row = characterize(events, p).iloc[0]
    assert row.during_coverage_fraction < 1
    p['valid_measurement'] = False
    p['value'] = np.nan
    row = characterize(events, p).iloc[0]
    assert row.status == 'insufficient_particle_coverage'
    assert np.isnan(row.during_vs_pre_pct)
    try:
        validate_particles(pd.concat([p, p.iloc[[0]]], ignore_index=True))
    except ValueError as exc:
        assert 'Duplicate' in str(exc)
    else:
        raise AssertionError('Duplicate measurement accepted')
    print('PASS: event overlap weighting, coverage gaps, no false zero baseline, no duplicated measurements')


def _digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def verify_inputs(rates_csv, events_csv):
    '''Require reader integrity and previously verified FIELDS strong-event catalog.'''
    reader_manifest = Path(str(rates_csv) + '.manifest.json')
    if not reader_manifest.is_file():
        raise ValueError(f'Missing particle reader manifest: {reader_manifest}')
    reader = json.loads(reader_manifest.read_text(encoding='utf-8'))
    if (reader.get('csv_sha256') != _digest(rates_csv) or
        reader.get('analysis_grade') != 'survey_only'):
        raise ValueError('Particle rate file checksum/product label mismatch')
    if not events_csv.name.endswith('_strong_events.csv'):
        raise ValueError('FIELDS input must be the strong_events.csv catalog, not all 3,701 rotations')
    retention = events_csv.with_name(events_csv.name.replace('_strong_events.csv', '_strong_retention.json'))
    if not retention.is_file():
        raise ValueError(f'Missing FIELDS strong-event retention manifest: {retention}')
    fields = json.loads(retention.read_text(encoding='utf-8'))
    if fields.get('strong_events_sha256') != _digest(events_csv):
        raise ValueError('FIELDS strong events file differs from verified retention manifest')


def _survey_main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--rates-csv', type=Path,
                        help='CSV generated by particle_data_reader.py')
    parser.add_argument('--events-csv', type=Path,
                        help='FIELDS strong_events.csv, never the 3,701 preliminary catalog')
    parser.add_argument('--context-minutes', type=float, default=30)
    parser.add_argument('--min-coverage', type=float, default=0.7)
    parser.add_argument('--min-samples', type=int, default=3)
    parser.add_argument('--output', type=Path,
                        default=Path('data/isois/fields_event_particle_survey.csv'))
    args = parser.parse_args()
    if args.self_test:
        _survey_self_test()
        return
    if not args.rates_csv or not args.events_csv:
        parser.error('--rates-csv and --events-csv required')
    if not args.rates_csv.is_file() or not args.events_csv.is_file():
        parser.error('Input CSV not found')
    try:
        verify_inputs(args.rates_csv, args.events_csv)
        particle = pd.read_csv(args.rates_csv)
        magnetic = pd.read_csv(args.events_csv)
        result = characterize(magnetic, particle,
                              context_minutes=args.context_minutes,
                              min_coverage=args.min_coverage,
                              min_samples=args.min_samples)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(args.output, index=False)
        manifest = {
            'fields_events_csv': str(args.events_csv), 'particle_rates_csv': str(args.rates_csv),
            'output': str(args.output), 'analysis_grade': 'survey_only',
            'context_minutes': args.context_minutes,
            'min_coverage': args.min_coverage, 'min_samples': args.min_samples,
            'rows': len(result),
            'status_counts': {str(k): int(v) for k, v in result.status.value_counts().items()},
            'warning': 'Exploratory product. Summary count rates are NOT calibrated differential flux. No physical correlation claims.',
        }
        args.output.with_suffix(args.output.suffix + '.manifest.json').write_text(
            json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
        print(f'[SAVED] {args.output}: {len(result)} event-variable rows')
        print('[STATUS]', manifest['status_counts'])
        print('[SURVEY ONLY] No statistical association or physical flux claim computed')
    except (ValueError, RuntimeError, OSError, KeyError) as exc:
        parser.exit(1, f'[ERROR] {exc}\n')

# === NATIVE 1–10s COVERAGE (from particle_coverage_audit.py) ===
import argparse
import json
import re
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd


PATTERN = re.compile(
    r'^psp_isois-epihi_l2-([\w-]+)_(\d{8})_v(\d+)'
    r'(_events_[0-9a-f]+)?\.npz$', re.I
)
NS = 1_000_000_000


def _inventory(root, products):
    """Choose latest CDF version per product/day; prefer full-day over subset."""
    out = {p: {} for p in products}
    for product in products:
        for file in (root / product).glob('*.npz'):
            m = PATTERN.fullmatch(file.name)
            if m is None or m.group(1).lower() != product:
                continue
            day, version = m.group(2), int(m.group(3))
            # For the same version, a full-day partition covers any event subset.
            rank = (version, int(m.group(4) is None), file.stat().st_size)
            previous = out[product].get(day)
            if previous is None or rank > previous[0]:
                out[product][day] = (rank, file)
    return {p: [x[1] for _, x in sorted(entries.items())]
            for p, entries in out.items()}


def _load_intervals(paths, *, verify_sha=True):
    starts, ends, failed = [], [], []
    for path in paths:
        try:
            manifest_path = Path(str(path) + '.manifest.json')
            if not manifest_path.is_file():
                raise ValueError('No provenance manifest')
            manifest = json.loads(manifest_path.read_text())
            if verify_sha and manifest.get('output_sha256') != _sha256(path):
                raise ValueError('Manifest checksum mismatch')
            with np.load(path, allow_pickle=False) as data:
                s = np.asarray(data['integration_start_utc_ns'], dtype=np.int64)
                e = np.asarray(data['integration_end_utc_ns'], dtype=np.int64)
            if s.ndim != 1 or e.ndim != 1 or s.shape != e.shape:
                raise ValueError('Invalid integration bounds arrays')
            if np.any(e <= s):
                raise ValueError('Nonpositive integration duration')
            if len(s):
                starts.append(s)
                ends.append(e)
        except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
            failed.append({'file': str(path),
                           'error': f'{type(exc).__name__}: {exc}'})
    if not starts:
        return np.empty(0, np.int64), np.empty(0, np.int64), failed
    s, e = np.concatenate(starts), np.concatenate(ends)
    order = np.lexsort((e, s))
    s, e = s[order], e[order]
    unique = np.r_[True, (s[1:] != s[:-1]) | (e[1:] != e[:-1])]
    return s[unique], e[unique], failed


def _segment_coverage(starts, ends, prefix_max_end, left, right):
    a, b = _utc(left).value, _utc(right).value
    if b <= a:
        raise ValueError('Invalid nonpositive event segment')
    if not len(starts):
        return dict(fraction=0.0, seconds=0.0, overlapping=0,
                    complete=0, midpoint=0)
    lo = int(np.searchsorted(prefix_max_end, a, side='right'))
    hi = int(np.searchsorted(starts, b, side='left'))
    s, e = starts[lo:hi], ends[lo:hi]
    valid = (s < b) & (e > a)
    s, e = s[valid], e[valid]
    if not len(s):
        return dict(fraction=0.0, seconds=0.0, overlapping=0,
                    complete=0, midpoint=0)
    lefts = np.maximum(s, a)
    rights = np.minimum(e, b)
    # Union of intervals: avoid double-credit for overlapping integrations.
    covered, x0, x1 = 0, None, None
    for si, ei in zip(lefts, rights):
        si, ei = int(si), int(ei)
        if x0 is None:
            x0, x1 = si, ei
        elif si <= x1:
            x1 = max(x1, ei)
        else:
            covered += x1 - x0
            x0, x1 = si, ei
    covered += x1-x0
    return dict(
        fraction=min(1., covered / (b-a)),
        seconds=covered/NS,
        overlapping=len(s),
        complete=int(np.count_nonzero((s >= a) & (e <= b))),
        midpoint=int(np.count_nonzero(((s + (e-s)//2) >= a) &
                                      ((s + (e-s)//2) < b))),
    )


def _audit(events, source_by_product, verify_sha):
    rows, inventory_info, bad_partitions = [], [], []
    for product, paths in source_by_product.items():
        starts, ends, failed = _load_intervals(paths, verify_sha=verify_sha)
        bad_partitions.extend({'product': product, **entry} for entry in failed)
        prefix = np.maximum.accumulate(ends) if len(ends) else ends
        med_seconds = (float(np.median((ends-starts)/NS))
                       if len(starts) else None)
        inventory_info.append({
            'product': product, 'partition_files_found': len(paths),
            'corrupt_or_unverifiable_files': len(failed),
            'distinct_time_integrations': len(starts),
            'median_integration_seconds': med_seconds,
        })
        for w in events.itertuples(index=False):
            start, stop = _utc(w.event_start_utc), _utc(w.event_end_exclusive_utc)
            duration = (stop-start).total_seconds()
            segments = {
                'before': _segment_coverage(
                    starts, ends, prefix, w.window_start_utc, start),
                'during': _segment_coverage(starts, ends, prefix, start, stop),
                'after': _segment_coverage(
                    starts, ends, prefix, stop, w.window_end_exclusive_utc),
            }
            # "Covered" means temporal overlap only; this cannot validate a
            # particle measurement or establish counts within an exact segment.
            during = segments['during']
            fine = med_seconds is not None and med_seconds <= 10.
            full = all(segments[k]['fraction'] >= .90
                       for k in ('before', 'during', 'after'))
            rows.append({
                'event_id': w.event_id, 'product': product,
                'event_duration_seconds': duration,
                'median_integration_seconds': med_seconds,
                'before_fraction_unfiltered': segments['before']['fraction'],
                'during_fraction_unfiltered': during['fraction'],
                'after_fraction_unfiltered': segments['after']['fraction'],
                'complete_integrations_during': during['complete'],
                'midpoints_during': during['midpoint'],
                'overlapping_integrations_during': during['overlapping'],
                'during_90pct_covered': during['fraction'] >= .90,
                'full_window_90pct_covered': full,
                'fine_10s_window_temporal_coverage': bool(
                    fine and full and during['complete'] >= 3),
                'precise_10s_flanks_approx_1s': bool(
                    med_seconds is not None and med_seconds <= 2 and full
                    and segments['before']['complete'] >= 5
                    and segments['after']['complete'] >= 5
                    and during['complete'] >= 3),
            })
    return pd.DataFrame(rows), inventory_info, bad_partitions


def _fine_self_test():
    t0 = _utc('2018-11-06T00:00:00Z').value
    s = t0 + np.arange(-10, 43, dtype=np.int64)*NS
    e = s + NS
    windows = pd.DataFrame([{
        'event_id': 'test', 'event_start_utc': '2018-11-06T00:00:00Z',
        'event_end_exclusive_utc': '2018-11-06T00:00:30Z',
        'window_start_utc': '2018-11-05T23:59:50Z',
        'window_end_exclusive_utc': '2018-11-06T00:00:40Z',
    }])
    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        folder = root / 'second-rates'
        folder.mkdir()
        path = folder / 'psp_isois-epihi_l2-second-rates_20181106_v26.npz'
        np.savez_compressed(
            path, integration_start_utc_ns=s, integration_end_utc_ns=e)
        Path(str(path)+'.manifest.json').write_text(json.dumps(
            {'output_sha256': _sha256(path)}))
        inv = _inventory(root, ['second-rates'])
        audit, info, errors = _audit(windows, inv, verify_sha=True)
        assert len(audit) == 1 and not errors
        row = audit.iloc[0]
        assert row['before_fraction_unfiltered'] == 1
        assert row['during_fraction_unfiltered'] == 1
        assert row['after_fraction_unfiltered'] == 1
        assert row['complete_integrations_during'] == 30
        assert row['fine_10s_window_temporal_coverage']
        assert row['precise_10s_flanks_approx_1s']
        # An absent 2-second interval reduces the actual union coverage.
        m = ~((s >= t0) & (s < t0+2*NS))
        p = np.maximum.accumulate(e[m])
        seg = _segment_coverage(s[m], e[m], p,
                                '2018-11-06T00:00:00Z',
                                '2018-11-06T00:00:30Z')
        assert abs(seg['fraction'] - 28/30) < 1e-12
    print('PASS: overlapping integration union, missing-sample gap, '
          '1s flanks, manifests and 10s event boundaries')


def _fine_main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--self-test', action='store_true')
    p.add_argument('--events-csv', type=Path)
    p.add_argument('--input-dir', type=Path,
                   default=Path('data/isois/highres'))
    p.add_argument('--output', type=Path,
                   default=Path('data/isois/highres/coverage_10s.csv'))
    p.add_argument('--before-seconds', type=float, default=10)
    p.add_argument('--after-seconds', type=float, default=10)
    p.add_argument('--no-verify-sha', action='store_true')
    args = p.parse_args(argv)
    if args.self_test:
        _fine_self_test()
        return
    if not args.events_csv:
        p.error('Provide --events-csv')
    if not np.isfinite(args.before_seconds) or not np.isfinite(args.after_seconds) \
            or min(args.before_seconds, args.after_seconds) <= 0:
        p.error('Audit requires positive pre/post padding')
    events = validate_fields_events(args.events_csv)
    windows = build_event_windows(
        events, before_minutes=args.before_seconds/60,
        after_minutes=args.after_seconds/60)
    inventory = _inventory(args.input_dir, DEFAULT_PRODUCTS)
    table, sources, failed = _audit(
        windows, inventory, verify_sha=not args.no_verify_sha)
    _atomic_csv(table, args.output)
    per_product = []
    for product, group in table.groupby('product', sort=False):
        per_product.append({
            'product': product,
            'events': int(len(group)),
            'during_90pct_covered': int(group['during_90pct_covered'].sum()),
            'full_window_90pct_covered': int(group['full_window_90pct_covered'].sum()),
            'fine_10s_window_temporal_coverage':
                int(group['fine_10s_window_temporal_coverage'].sum()),
            'precise_10s_flanks_approx_1s':
                int(group['precise_10s_flanks_approx_1s'].sum()),
        })
    any_fine = table.groupby('event_id')['fine_10s_window_temporal_coverage'].any()
    any_precise = table.groupby('event_id')['precise_10s_flanks_approx_1s'].any()
    summary = {
        'evaluation': 'NOMINAL TIME COVERAGE ONLY; not particle-quality or counting-statistics adequacy',
        'events': len(events),
        'padding_seconds': {'before': args.before_seconds,
                            'after': args.after_seconds},
        'nominal_fine_window_covered_by_any_product': int(any_fine.sum()),
        'approx_1s_flanks_covered_by_any_product': int(any_precise.sum()),
        'product_summaries': per_product,
        'source_inventory': sources,
        'invalid_or_corrupt_partitions': failed,
        'quality_flags_interpreted': False,
        'particle_counts_or_poisson_significance_checked': False,
    }
    summary_path = args.output.with_suffix('.summary.json')
    _atomic_text(summary_path, json.dumps(summary, indent=2) + '\n')
    print(f'[COVERAGE] {len(events)} FIELDS events; '
          f'{args.before_seconds:g}s before and {args.after_seconds:g}s after')
    for item in per_product:
        n = item['events']
        print(f'  {item["product"]}: during >=90% {item["during_90pct_covered"]}/{n}; '
              f'whole window >=90% {item["full_window_90pct_covered"]}/{n}; '
              f'fine-cadence window {item["fine_10s_window_temporal_coverage"]}/{n}')
    print('[TOTAL] Nominal 1-10s full-window coverage for '
          f'{any_fine.sum()}/{len(events)} candidates')
    print('[TOTAL] Approximate 1s precise-flank coverage for '
          f'{any_precise.sum()}/{len(events)} candidates')
    print(f'[SAVED] {args.output} and {summary_path}')
    print('[CAUTION] Quality flags, energy-specific availability, '
          'particle-count statistics and causal association remain UNTESTED.')

# === PER-CHANNEL MULTIRESOLUTION COVERAGE (from particle_multires_coverage.py) ===
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


# Strict same telescope, side and particle species for direct fallback.
# "rate" / "flux" and energy ranges STILL differ between products.
CHANNELS = {
    'LET1_A_H': {
        'second-rates': ('LET1_A_H_Rate',),
        'let1-rates10': ('A_H_Flux', 'A_H_Rate'),
        'let1-rates60': ('A_H_Flux', 'A_H_Rate'),
        'let1-rates3600': ('A_H_Flux', 'A_H_Rate'),
        'summary': ('A_H_Rate_TS',),
    },
    'HET_A_H': {
        'second-rates': ('HET_A_H_Rate',),
        'het-rates10': ('A_H_Flux', 'A_H_Rate'),
        'het-rates60': ('A_H_Flux', 'A_H_Rate'),
        'het-rates3600': ('A_H_Flux', 'A_H_Rate'),
        'summary': ('HET_A_H_Rate_TS',),
    },
    'HET_A_Electrons': {
        'second-rates': ('HET_A_Electrons_Rate',),
        'het-rates10': ('A_Electrons_Flux', 'A_Electrons_Rate'),
        'het-rates60': ('A_Electrons_Flux', 'A_Electrons_Rate'),
        'het-rates3600': ('A_Electrons_Flux', 'A_Electrons_Rate'),
        'summary': ('HET_A_Electrons_Rate_TS',),
    },
    # The combined Z=6-28 rate has NO exact species-by-species equivalent.
    'LET1_A_Heavy_Z6_28': {
        'summary': ('A_Heavy_Rate_TS',),
    },
    # Five-minute L2 sectored rates provide distinct HEAVY-ION GROUPS.
    # They are NOT the summary's combined Z=6..28 rate; preserve their
    # telescope, range (R1/R26), side (A) and grouping identity.
    'LET1_R1A_CNO': {'let1-rates300': ('R1A_CNO_SECT_Rate',)},
    'LET1_R1A_NetoSi': {'let1-rates300': ('R1A_NetoSi_SECT_Rate',)},
    'LET1_R1A_FeGroup': {'let1-rates300': ('R1A_FeGroup_SECT_Rate',)},
    'LET1_R26A_CNO': {'let1-rates300': ('R26A_CNO_SECT_Rate',)},
    'LET1_R26A_NetoSi': {'let1-rates300': ('R26A_NetoSi_SECT_Rate',)},
    'LET1_R26A_FeGroup': {'let1-rates300': ('R26A_FeGroup_SECT_Rate',)},
    'HET_A_CNO': {'het-rates300': ('A_CNO_SECT_Rate',)},
    'HET_A_NetoSi': {'het-rates300': ('A_NetoSi_SECT_Rate',)},
    'HET_A_FeGroup': {'het-rates300': ('A_FeGroup_SECT_Rate',)},
}
HIGHR_PRODUCTS = tuple(sorted({
    product for options in CHANNELS.values()
    for product in options if product != 'summary'
}))
PRODUCT_CADENCE = {
    'second-rates': 1,
    'let1-rates10': 10, 'het-rates10': 10,
    'let1-rates60': 60, 'let2-rates60': 60, 'het-rates60': 60,
    'let1-rates300': 300, 'het-rates300': 300,
    'let1-rates3600': 3600, 'het-rates3600': 3600,
    'summary': 3600,
}
TIER_PRIORITY = {
    'fine_event_and_flanks': 0,
    'fine_event_only': 1,
    'coarse_event_scale': 2,
    'fine_partial_context': 3,
    'coarse_partial_context': 4,
    'hourly_background_only': 5,
    'no_observation': 6,
}


def _select_paths_with_channel(paths, names):
    """Only consider files whose own manifest actually contains this channel."""
    matching, rejected = [], []
    for path in paths:
        m = Path(str(path) + '.manifest.json')
        try:
            manifest = json.loads(m.read_text())
            variables = manifest['variables']
            if not isinstance(variables, dict):
                raise ValueError('variables metadata is not a mapping')
            found = next((name for name in names if name in variables), None)
            if found:
                matching.append((path, found))
        except (OSError, ValueError, KeyError, TypeError) as exc:
            rejected.append({'file': str(path), 'problem': str(exc)})
    return matching, rejected


def _summary_sources(path):
    """Use hourly *survey* rates as coarse data; never call them calibrated flux."""
    by_channel = {}
    if path is None or not path.is_file():
        return by_channel, None
    frame = pd.read_csv(path)
    required = {
        'variable', 'integration_start_utc', 'integration_end_utc',
        'valid_measurement'
    }
    if required - set(frame.columns):
        raise ValueError(f'Summary CSV missing {required - set(frame.columns)}')
    valid = frame['valid_measurement'].astype(str).str.lower().isin(
        ('true', '1', 'yes'))
    frame = frame.loc[valid].copy()
    for channel, product_vars in CHANNELS.items():
        name = product_vars.get('summary', (None,))[0]
        group = frame.loc[frame['variable'] == name]
        if group.empty:
            continue
        # Pandas may parse ISO timestamps at microsecond resolution on some
        # versions. Normalizing to ns is essential: FIELDS event bounds are ns.
        # Without this, hourly measurements appear to have zero overlap.
        start = pd.DatetimeIndex(
            pd.to_datetime(group['integration_start_utc'], utc=True, errors='raise')
        ).as_unit('ns').asi8.copy()
        end = pd.DatetimeIndex(
            pd.to_datetime(group['integration_end_utc'], utc=True, errors='raise')
        ).as_unit('ns').asi8.copy()
        if np.any(end <= start):
            raise ValueError(f'{name}: summary integration end <= start')
        indices = np.lexsort((end, start))
        s, e = start[indices], end[indices]
        distinct = np.r_[True, (s[1:] != s[:-1]) | (e[1:] != e[:-1])]
        by_channel[channel] = (s[distinct], e[distinct], name)
    return by_channel, _sha256(path)


def _measure(event, starts, ends, *, source, variable, expected_cadence):
    begin, finish = _utc(event.event_start_utc), _utc(event.event_end_exclusive_utc)
    prefix = np.maximum.accumulate(ends) if len(ends) else ends
    pre = _segment_coverage(
        starts, ends, prefix, event.window_start_utc, begin)
    during = _segment_coverage(starts, ends, prefix, begin, finish)
    post = _segment_coverage(
        starts, ends, prefix, finish, event.window_end_exclusive_utc)
    # Use the actual stored integration width, not filename cadence, when possible.
    widths = (ends-starts)/NS
    cadence = (float(np.median(widths)) if len(widths)
               else float(expected_cadence))
    nominal_full = all(x['fraction'] >= .90 for x in (pre, during, post))
    if (cadence <= 10 and nominal_full
            and during['complete'] >= 3):
        tier = 'fine_event_and_flanks'
    elif cadence <= 10 and during['fraction'] >= .90:
        tier = 'fine_event_only'
    elif cadence <= 10 and during['fraction'] > 0:
        tier = 'fine_partial_context'
    elif (cadence <= 300 and during['fraction'] >= .90
          and during['complete'] >= 1):
        tier = 'coarse_event_scale'
    elif cadence <= 300 and during['fraction'] > 0:
        tier = 'coarse_partial_context'
    elif during['fraction'] > 0:
        tier = 'hourly_background_only'
    else:
        tier = 'no_observation'
    clean = (
        cadence <= 10 and during['complete'] >= 3
        and pre['complete'] >= 1 and post['complete'] >= 1
        and nominal_full
    )
    return {
        'event_id': event.event_id,
        'product': source,
        'measurement_variable': variable,
        'integration_seconds': cadence,
        'before_fraction_unfiltered': pre['fraction'],
        'during_fraction_unfiltered': during['fraction'],
        'after_fraction_unfiltered': post['fraction'],
        'full_window_90pct_nominal': nominal_full,
        'complete_integrations_during': during['complete'],
        'complete_integrations_before': pre['complete'],
        'complete_integrations_after': post['complete'],
        'isolated_before_after_integrations': bool(clean),
        'coverage_tier': tier,
        'is_calibrated_flux': variable.endswith('_Flux'),
        'measurement_units_type': (
            'differential_flux' if variable.endswith('_Flux')
            else 'counts_per_second'),
        'independent_of_other_products': False,
        'instrument_quality_validated': False,
        'counts_significance_validated': False,
    }


def evaluate(events, inventory, summary_sources, *, verify_sha=True):
    """Return per-source rows and source inventory; never sum measurements."""
    rows, source_rows, errors = [], [], []
    # Load and decode time intervals only once per product+channel.
    for channel, sources in CHANNELS.items():
        for product, candidates in sources.items():
            if product == 'summary':
                if channel not in summary_sources:
                    continue
                start, end, actual_name = summary_sources[channel]
                matched_files = []
                bad = []
            else:
                matched_files, bad = _select_paths_with_channel(
                    inventory.get(product, []), candidates)
                selected = [p for p, _ in matched_files]
                start, end, corrupt = _load_intervals(
                    selected, verify_sha=verify_sha)
                bad += [
                    {'file': item['file'], 'problem': item['error']}
                    for item in corrupt
                ]
                actual_name = (
                    matched_files[0][1] if matched_files else
                    candidates[0]
                )
            errors.extend(
                {'channel': channel, 'product': product, **e} for e in bad
            )
            if not len(start):
                continue
            source_rows.append({
                'channel': channel, 'product': product,
                'actual_variable': actual_name,
                'partitions': len(matched_files) if product != 'summary' else None,
                'integration_records': int(len(start)),
                'median_integration_seconds':
                    float(np.median((end-start)/NS)),
                'calibrated_flux': actual_name.endswith('_Flux'),
                'note': ('Summary is count-rate survey, not flux' if product == 'summary'
                         else 'Each product has its own energy bins and quality flags'),
            })
            for ev in events.itertuples(index=False):
                record = _measure(ev, start, end, source=product,
                                  variable=actual_name,
                                  expected_cadence=PRODUCT_CADENCE[product])
                if record['coverage_tier'] != 'no_observation':
                    rows.append({'channel': channel, **record})
    return pd.DataFrame(rows), pd.DataFrame(source_rows), errors


def choose_best(events, records):
    """Single preferred provenance per event and channel, not a summed signal."""
    base = pd.MultiIndex.from_product(
        [[e.event_id for e in events.itertuples(index=False)], CHANNELS],
        names=['event_id', 'channel']).to_frame(index=False)
    if records.empty:
        base['best_tier'] = 'no_observation'
        base['best_product'] = ''
        base['best_variable'] = ''
        base['nominal_full_window'] = False
        base['alternative_sources_available'] = 0
        return base
    r = records.copy()
    r['_tier_rank'] = r['coverage_tier'].map(TIER_PRIORITY)
    r['_source_rank'] = r['integration_seconds']
    # Sorted tier then cadence; at equal cadence prefer a calibrated flux
    # (but these remain different energy bins and are not equivalent).
    r = r.sort_values(
        ['event_id', 'channel', '_tier_rank', '_source_rank', 'is_calibrated_flux'],
        ascending=[True, True, True, True, False], kind='stable')
    best = r.drop_duplicates(['event_id', 'channel'])
    counts = (r.groupby(['event_id', 'channel'], sort=False).size()
              .rename('alternative_sources_available').reset_index())
    output = base.merge(best[[
        'event_id', 'channel', 'coverage_tier', 'product',
        'measurement_variable', 'full_window_90pct_nominal',
        'integration_seconds', 'isolated_before_after_integrations'
    ]], on=['event_id', 'channel'], how='left')
    output = output.merge(counts, on=['event_id', 'channel'], how='left')
    output = output.rename(columns={
        'coverage_tier': 'best_tier',
        'product': 'best_product',
        'measurement_variable': 'best_variable',
        'full_window_90pct_nominal': 'nominal_full_window',
        'integration_seconds': 'best_integration_seconds',
    })
    output['best_tier'] = output['best_tier'].fillna('no_observation')
    output['best_product'] = output['best_product'].fillna('')
    output['best_variable'] = output['best_variable'].fillna('')
    output['nominal_full_window'] = output['nominal_full_window'].fillna(False).astype(bool)
    output['alternative_sources_available'] = (
        output['alternative_sources_available'].fillna(0).astype(int))
    return output


def _multires_self_test():
    from tempfile import TemporaryDirectory

    windows = pd.DataFrame([
        {'event_id': 'a', 'event_start_utc': '2018-11-06T00:00:00Z',
         'event_end_exclusive_utc': '2018-11-06T00:03:00Z',
         'window_start_utc': '2018-11-05T23:59:50Z',
         'window_end_exclusive_utc': '2018-11-06T00:03:10Z'},
        {'event_id': 'b', 'event_start_utc': '2018-11-06T01:00:00Z',
         'event_end_exclusive_utc': '2018-11-06T01:02:00Z',
         'window_start_utc': '2018-11-06T00:59:50Z',
         'window_end_exclusive_utc': '2018-11-06T01:02:10Z'}
    ])
    start = _utc('2018-11-06T00:00:00Z').value
    # 1-second before/during/after with all but 1 sample.
    s = start + np.arange(-10, 191, dtype=np.int64)*NS
    e = s+NS
    fine = _measure(windows.iloc[0], s, e, source='second-rates',
                    variable='LET1_A_H_Rate', expected_cadence=1)
    assert fine['coverage_tier'] == 'fine_event_and_flanks'
    assert fine['isolated_before_after_integrations']
    # Fine cadence with partial flanks is still fine data, not coarser context.
    shortened = _measure(
        windows.iloc[0], s[(s >= start) & (s < start+180*NS)],
        e[(s >= start) & (s < start+180*NS)],
        source='second-rates', variable='LET1_A_H_Rate',
        expected_cadence=1)
    assert shortened['coverage_tier'] == 'fine_event_only'
    # Hourly integration covers the whole time range, but NEVER the 10s flanks.
    hourly = _measure(
        windows.iloc[0],
        np.asarray([start-1800*NS]),
        np.asarray([start+1800*NS]),
        source='summary', variable='A_H_Rate_TS', expected_cadence=3600)
    assert hourly['coverage_tier'] == 'hourly_background_only'
    assert not hourly['isolated_before_after_integrations']
    # One 60 s integration spanning event edge is context only.
    coarse = _measure(
        windows.iloc[0],
        np.asarray([start-30*NS]),
        np.asarray([start+30*NS]),
        source='let1-rates60', variable='A_H_Flux',
        expected_cadence=60)
    assert coarse['coverage_tier'] == 'coarse_partial_context'
    with TemporaryDirectory() as folder:
        survey = Path(folder) / 'summary.csv'
        pd.DataFrame([{
            'variable': 'A_H_Rate_TS',
            'integration_start_utc': '2018-11-05T23:30:00Z',
            'integration_end_utc': '2018-11-06T00:30:00Z',
            'valid_measurement': True,
        }, {
            'variable': 'HET_A_H_Rate_TS',
            'integration_start_utc': '2018-11-05T23:30:00Z',
            'integration_end_utc': '2018-11-06T00:30:00Z',
            'valid_measurement': False,
        }]).to_csv(survey, index=False)
        summary, sh = _summary_sources(survey)
        assert summary['LET1_A_H'][0][0] == _utc('2018-11-05T23:30:00Z').value
        assert summary['LET1_A_H'][1][0] == _utc('2018-11-06T00:30:00Z').value
        assert sh == _sha256(survey)
        assert 'LET1_A_H' in summary
        assert 'HET_A_H' not in summary
        fallback = _measure(windows.iloc[0], summary['LET1_A_H'][0],
                            summary['LET1_A_H'][1], source='summary',
                            variable='A_H_Rate_TS', expected_cadence=3600)
        assert fallback['coverage_tier'] == 'hourly_background_only'
    # A 300s CNO observation is a separate channel and must never
    # replace a proton spectrum or the combined summary Z6-28 group.
    with TemporaryDirectory() as folder:
        root = Path(folder)
        p = root / 'let1-rates300'
        p.mkdir()
        path = p/'psp_isois-epihi_l2-let1-rates300_20181106_v26.npz'
        np.savez_compressed(
            path,
            integration_start_utc_ns=np.asarray([start-50*NS]),
            integration_end_utc_ns=np.asarray([start+250*NS]),
            R26A_CNO_SECT_Rate=np.ones((1,2,25)))
        Path(str(path)+'.manifest.json').write_text(
            json.dumps({'output_sha256': _sha256(path), 'variables': {
                'R26A_CNO_SECT_Rate': {'units':'s^-1'}
            }}))
        inv = _inventory(root, HIGHR_PRODUCTS)
        rates, sources, failed = evaluate(windows, inv, {}, verify_sha=True)
        assert not failed
        assert len(sources) == 1
        assert set(rates['channel']) == {'LET1_R26A_CNO'}
        assert rates.iloc[0]['product'] == 'let1-rates300'
        assert rates.iloc[0]['measurement_units_type'] == 'counts_per_second'
        assert rates.iloc[0]['coverage_tier'] == 'coarse_partial_context'
        best = choose_best(windows, rates)
        assert best[(best.event_id=='a') &
                    (best.channel=='LET1_A_Heavy_Z6_28')].iloc[0]['best_tier'] == 'no_observation'
        assert best[(best.event_id=='a') &
                    (best.channel=='LET1_A_H')].iloc[0]['best_tier'] == 'no_observation'
    with TemporaryDirectory() as folder:
        root = Path(folder)
        p = root/'let1-rates3600'
        p.mkdir()
        path = p/'psp_isois-epihi_l2-let1-rates3600_20181106_v26.npz'
        np.savez_compressed(path, integration_start_utc_ns=np.asarray([start-1800*NS]),
                            integration_end_utc_ns=np.asarray([start+1800*NS]))
        Path(str(path)+'.manifest.json').write_text(
            json.dumps({'output_sha256':_sha256(path), 'variables':{
                'A_H_Flux':{'units':'flux'}, 'A_He_Flux':{'units':'flux'}
            }}))
        inv = _inventory(root, HIGHR_PRODUCTS)
        records, source_records, failures = evaluate(
            windows, inv, {}, verify_sha=True)
        assert not failures and len(source_records) == 1
        assert (records['channel'] == 'LET1_A_H').all()
        assert records.iloc[0]['coverage_tier'] == 'hourly_background_only'
        best = choose_best(windows, records)
        assert len(best) == len(windows)*len(CHANNELS)
        assert best[(best.event_id=='b')&(best.channel=='LET1_A_H')].iloc[0]['best_tier'] == 'no_observation'
        # The same detector with 1s data is preferred, never added to 1h.
        sampled = pd.concat([
            records,
            pd.DataFrame([{'channel':'LET1_A_H', **fine}])
        ], ignore_index=True)
        chosen = choose_best(windows, sampled)
        row = chosen[(chosen.event_id=='a')&(chosen.channel=='LET1_A_H')].iloc[0]
        assert row['best_product'] == 'second-rates'
        assert row['alternative_sources_available'] == 2
    print('PASS: exact channel identity, 1s flanks, 60s event context, '
          'hourly background, manifests and no double counting')


def _multires_main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--self-test', action='store_true')
    p.add_argument('--events-csv', type=Path)
    p.add_argument('--input-dir', type=Path,
                   default=Path('data/isois/highres'))
    p.add_argument('--summary-csv', type=Path,
                   default=Path('data/isois/fields_event_summary_rates.csv'))
    p.add_argument('--output-dir', type=Path,
                   default=Path('data/isois/highres'))
    p.add_argument('--before-seconds', type=float, default=10)
    p.add_argument('--after-seconds', type=float, default=10)
    p.add_argument('--no-verify-sha', action='store_true')
    args = p.parse_args(argv)
    if args.self_test:
        _multires_self_test()
        return
    if not args.events_csv:
        p.error('--events-csv is required')
    if (not np.isfinite(args.before_seconds) or
            not np.isfinite(args.after_seconds) or
            min(args.before_seconds, args.after_seconds) <= 0):
        p.error('Require positive before/after padding')
    events = validate_fields_events(args.events_csv)
    windows = build_event_windows(
        events, before_minutes=args.before_seconds/60,
        after_minutes=args.after_seconds/60)
    inv = _inventory(args.input_dir, HIGHR_PRODUCTS)
    summary, summary_sha = _summary_sources(args.summary_csv)
    observations, sources, bad = evaluate(
        windows, inv, summary, verify_sha=not args.no_verify_sha)
    best = choose_best(windows, observations)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    _atomic_csv(best, args.output_dir/'multires_best_per_event_channel.csv')
    _atomic_csv(observations, args.output_dir/'multires_observed_sources.csv')
    _atomic_csv(sources, args.output_dir/'multires_source_inventory.csv')
    counts = (best.groupby(['channel','best_tier']).size().unstack(fill_value=0)
              .reindex(columns=list(TIER_PRIORITY), fill_value=0))
    per_channel = counts.reset_index().to_dict(orient='records')
    # All four original channels are independent; show totals PER channel,
    # never sum their observations across species/detectors.
    report = {
        'fields_candidates': len(events),
        'pad_seconds': {'before':args.before_seconds, 'after':args.after_seconds},
        'interpretation': (
            'COVERAGE TIMING ONLY. Flux and count rates are NOT merged. '
            'No reconstruction of missing 1s/10s measurements. '
            'No quality flags, energy compatibility or Poisson significance verified.'
        ),
        'per_channel_coverage_tiers': per_channel,
        'source_inventory': sources.to_dict('records'),
        'rejected_source_partitions': bad,
        'summary_survey': str(args.summary_csv),
        'summary_survey_sha256': summary_sha,
        'metadata_quality_flags_interpreted': False,
        'channel_energy_spectra_verified_compatible': False,
        'particle_counts_statistically_adequate': None,
        'heavy_ion_warning': (
            'LET1 summary heavy Z6-28 rate and 5-minute LET1/HET sectored '
            'CNO/NetoSi/FeGroup count rates have DIFFERENT species groups, '
            'energy selections and viewing geometry. Never merge them as '
            'equivalent calibrated differential flux.'
        ),
    }
    output = args.output_dir/'multires_coverage_summary.json'
    _atomic_text(output, json.dumps(report, indent=2)+'\n')
    print(f'[MULTIRES] {len(events)} FIELDS candidates, '
          f'{args.before_seconds:g}s before/after = '
          f'{args.before_seconds:g}/{args.after_seconds:g}s')
    for row in per_channel:
        print('  '+row['channel']+': '
              +', '.join(k+'='+str(row[k]) for k in TIER_PRIORITY))
    print(f'[FILES] {len(sources)} exact-channel source groups; '
          f'{len(bad)} invalid/rejected partitions')
    print(f'[SAVED] {output} and multires_*.csv')
    print('[CAUTION] Coverage tiers are TIMING/CONTEXT only. '
          'Do not equate 3600s background with 10s event resolution.')


# ---------------------------------------------------------------------------
# Unified public analysis CLI (all algorithms are implemented above)
# ---------------------------------------------------------------------------
def _default_coverage_args(mode, args):
    from pathlib import Path
    args = list(args)
    if any(x in args for x in ("--self-test", "--help", "-h")):
        return args
    def has_flag(flag):
        return any(arg == flag or arg.startswith(flag + "=") for arg in args)
    if mode in ("coverage", "fine-coverage") and not has_flag("--events-csv"):
        path = Path("data/full_analysis_nov_dec_2018_candidates_strong_events.csv")
        if path.is_file():
            args.extend(["--events-csv", str(path)])
    if mode == "coverage" and not has_flag("--summary-csv"):
        path = Path("data/isois/fields_event_summary_rates.csv")
        if path.is_file():
            args.extend(["--summary-csv", str(path)])
    return args


def self_test():
    _survey_self_test()
    _fine_self_test()
    _multires_self_test()
    assert _default_coverage_args("coverage", ["--self-test"]) == ["--self-test"]
    assert _default_coverage_args("coverage", ["--events-csv=custom.csv"])[:1] == [
        "--events-csv=custom.csv"]
    print("PASS: survey, fine cadence, multiresolution analysis and unified CLI")


def help_text():
    print("""IS⊙IS analysis — all analysis logic lives in this one Python file

  python project.py isois-analyzer coverage [options]
  python project.py isois-analyzer fine-coverage [options]
  python project.py isois-analyzer survey --rates-csv FILE --events-csv FILE
  python project.py isois-analyzer --self-test

Append --help to any mode to view its full input/output settings.
Existing results still use the exact same CSV/JSON formats and paths.
Quality/Poisson significance remain unverified; hourly is not fine timing.
""")


def main(argv=None):
    import sys
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("--help", "-h", "help"):
        help_text()
        return
    if argv == ["--self-test"]:
        self_test()
        return
    mode, *args = argv
    if mode == "coverage":
        return _multires_main(_default_coverage_args(mode, args))
    if mode == "fine-coverage":
        return _fine_main(_default_coverage_args(mode, args))
    if mode == "survey":
        original = sys.argv
        try:
            sys.argv = [original[0], *args]
            return _survey_main()
        finally:
            sys.argv = original
    raise SystemExit(f"Unknown analysis mode {mode!r}; choose coverage, fine-coverage or survey")


if __name__ == "__main__":
    main()
