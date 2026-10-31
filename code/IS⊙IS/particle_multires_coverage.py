#!/usr/bin/env python3
"""Same-channel multiresolution EPI-Hi coverage audit for FIELDS switchbacks.

Keeps detector/species/view identity. This program DOES NOT interpolate,
combine counts from overlapping products, or claim that different energy
bins have the same physical response. Quality flags and Poisson significance
remain separate scientific checks.

Three useful tiers: 1-10 s fine timing (including +/-10 s flanks);
60-300 s event-scale context; 3600 s hourly background only.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from particle_data_reader import (
    _atomic_csv, _atomic_text, _sha256, _utc,
    build_event_windows, validate_fields_events,
)
from particle_coverage_audit import (
    _inventory, _load_intervals, _segment_coverage, NS,
)

# Strict same telescope, side and particle species for direct fallback.
# "rate" / "flux" and energy ranges STILL differ between products.
CHANNELS = {
    'LET1_A_H': {
        'second-rates': ('LET1_A_H_Rate',),
        'let1-rates10': ('A_H_Flux', 'A_H_Rate'),
        'let1-rates60': ('A_H_Flux', 'A_H_Rate'),
        'let1-rates300': ('A_H_Flux', 'A_H_Rate'),
        'let1-rates3600': ('A_H_Flux', 'A_H_Rate'),
        'summary': ('A_H_Rate_TS',),
    },
    'HET_A_H': {
        'second-rates': ('HET_A_H_Rate',),
        'het-rates10': ('A_H_Flux', 'A_H_Rate'),
        'het-rates60': ('A_H_Flux', 'A_H_Rate'),
        'het-rates300': ('A_H_Flux', 'A_H_Rate'),
        'het-rates3600': ('A_H_Flux', 'A_H_Rate'),
        'summary': ('HET_A_H_Rate_TS',),
    },
    'HET_A_Electrons': {
        'second-rates': ('HET_A_Electrons_Rate',),
        'het-rates10': ('A_Electrons_Flux', 'A_Electrons_Rate'),
        'het-rates60': ('A_Electrons_Flux', 'A_Electrons_Rate'),
        'het-rates300': ('A_Electrons_Flux', 'A_Electrons_Rate'),
        'het-rates3600': ('A_Electrons_Flux', 'A_Electrons_Rate'),
        'summary': ('HET_A_Electrons_Rate_TS',),
    },
    # The combined Z=6-28 rate has NO exact species-by-species equivalent.
    'LET1_A_Heavy_Z6_28': {
        'summary': ('A_Heavy_Rate_TS',),
    },
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


def self_test():
    from tempfile import TemporaryDirectory
    from particle_data_reader import _sha256

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


def main(argv=None):
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
        self_test()
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
            'LET1 heavy Z6-28 summary group cannot be replaced as same '
            'channel by separate LET1 elemental or grouped heavy-ion flux.'
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


if __name__ == '__main__':
    main()
