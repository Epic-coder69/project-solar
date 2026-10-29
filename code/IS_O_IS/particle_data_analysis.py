#!/usr/bin/env python3
"""ISʘIS survey-rate characterization centered on FIELDS magnetic events.

This is intentionally *not* the final field-particle correlation engine.
Measures independent pre/during/post count-rate summaries from the L2 summary
survey product. No p-values, flux interpretation, energy spectra or causal claims.
All windows use actual exposure overlap, not nearest-minute timestamp matching.
"""
from __future__ import annotations

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


def _utc(series):
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
        frame[col] = _utc(frame[col])
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


def self_test():
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


def main():
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
        self_test()
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


if __name__ == '__main__':
    main()
