#!/usr/bin/env python3
"""Audit true-time EPI-Hi NPZ coverage of FIELDS switchback candidates.

IMPORTANT: This is a NOMINAL TIMING audit, NOT a scientific validity audit.
It does not interpret Quality_Flag, Poisson uncertainties, counts, or flux.
Different cadence products can share the underlying detected particles; this
script never sums counts/rates/flux or invents missing fine-resolution samples.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from particle_data_reader import (
    _atomic_csv, _atomic_text, _sha256, _utc,
    build_event_windows, validate_fields_events,
)
from particle_highres_reader import DEFAULT_PRODUCTS

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


def self_test():
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


def main(argv=None):
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
        self_test()
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


if __name__ == '__main__':
    main()
