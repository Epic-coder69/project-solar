"""Unattended, storage-aware FIELDS switchback pipeline.

Run the reader and magnetic analysis one UTC day at a time. Before retiring
any full-cadence CDF, archive ALL detected events at original sample cadence
and all analyzed one-second candidate windows with their computed baselines.

Deletion is opt-in (--cleanup-raw), limited to canonical FIELDS cache files
inside data/fields_full whose days have completed verified archives.
Original NASA CDFs can be fetched again if scientifically needed.
"""
import argparse
import gzip
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from magnetic_data_reader import (FULL_CACHE, FULL_NAME, _parse_duration,
                                  _as_utc_naive, read_full_interval)
from magnetic_data_analysis import (FULL_EVENT_COLUMNS,
                                    EVENT_REGISTRY_COLUMNS, deduplicate)

DATA = Path('data')
COLUMNS = ['BR', 'BT', 'BN', 'Bmag']


def atomic_csv(frame, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.tmp')
    try:
        frame.to_csv(temp, index=False)
        temp.replace(path)
    finally:
        temp.unlink(missing_ok=True)


def atomic_json(data, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.tmp')
    try:
        temp.write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')
        temp.replace(path)
    finally:
        temp.unlink(missing_ok=True)


def sha256_file(path):
    checksum = hashlib.sha256()
    with Path(path).open('rb') as handle:
        while True:
            block = handle.read(1024 * 1024)
            if not block:
                break
            checksum.update(block)
    return checksum.hexdigest()


def validate_gzip(path):
    """Read entire compressed stream to check its CRC and trailer (bounded RAM)."""
    with gzip.open(path, 'rb') as handle:
        while handle.read(1024 * 1024):
            pass


def archive_native_event(event, archive_dir, context_minutes=5,
                         stream_reader=read_full_interval):
    """Persist original-cadence, quality-filtered RTN vectors in gzip CSV.

    Includes five minutes of surrounding data. This preserves measured vector
    samples, not a byte-for-byte copy of NASA CDF instrument metadata/flags.
    """
    start = pd.Timestamp(event['event_start_utc'])
    end_inclusive = pd.Timestamp(event['event_end_utc'])
    if end_inclusive < start:
        raise ValueError('Event boundaries out of order')
    if int(event.get('native_valid_samples', 0)) <= 0:
        raise ValueError('Event has no validated native magnetic samples')
    left = start - pd.Timedelta(minutes=context_minutes)
    right = end_inclusive + pd.Timedelta(seconds=1, minutes=context_minutes)
    archive_dir = Path(archive_dir)
    archive_dir.mkdir(parents=True, exist_ok=True)
    name = str(event['event_id'])
    if not name.startswith('PSP-MAG-') or not all(
            ch.isalnum() or ch in '-_' for ch in name):
        raise ValueError('Unsafe event ID: ' + name)
    path = archive_dir / f'{name}.csv.gz'
    temporary = path.with_name(path.name + '.part')
    count = 0
    inside_event = 0
    valid_inside_event = 0
    first = True
    event_end_exclusive = end_inclusive + pd.Timedelta(seconds=1)
    try:
        with gzip.open(temporary, 'wt', encoding='utf-8', newline='') as out:
            for chunk in stream_reader(left, right):
                chunk = chunk.loc[(chunk.index >= left) &
                                  (chunk.index < right), COLUMNS]
                if chunk.empty:
                    continue
                in_event = chunk.loc[(chunk.index >= start) &
                                     (chunk.index < event_end_exclusive)]
                inside_event += len(in_event)
                valid_inside_event += int(
                    np.isfinite(in_event[['BR','BT','BN']]
                                .to_numpy(dtype=float)).all(axis=1).sum())
                chunk.to_csv(out, header=first, index_label='UTC')
                count += len(chunk)
                first = False
        if (first or inside_event == 0 or
                valid_inside_event < int(event['native_valid_samples'])):
            raise RuntimeError(f'Archived event incomplete for {name}: '
                               f'{valid_inside_event} valid samples in event, '
                               f'expected {event["native_valid_samples"]}')
        validate_gzip(temporary)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
    return dict(event_id=name, archive_file=str(path),
                first_requested_utc=left.isoformat(),
                end_requested_exclusive_utc=right.isoformat(),
                native_archived_samples=count,
                native_event_archived_samples=inside_event,
                native_event_valid_samples=valid_inside_event,
                sha256=sha256_file(path))


def analysis_outputs(catalog):
    name = catalog.stem
    prefix = DATA / f'full_analysis_{name}'
    return (Path(str(prefix) + '_events.csv'),
            Path(str(prefix) + '_registry.csv'),
            Path(str(prefix) + '_audit.csv'))


def validate_analysis(candidate_csv, events_path, registry_path, audit_path):
    """Do not approve archival/deletion when ANY candidate failed processing."""
    seeds = pd.read_csv(candidate_csv)
    events = pd.read_csv(events_path) if not seeds.empty else pd.DataFrame(
        columns=FULL_EVENT_COLUMNS)
    registry = pd.read_csv(registry_path) if not seeds.empty else pd.DataFrame(
        columns=EVENT_REGISTRY_COLUMNS)
    audit = pd.read_csv(audit_path) if not seeds.empty else pd.DataFrame()
    if not seeds.empty:
        if len(audit) != len(seeds):
            raise RuntimeError('Audit does not cover all candidate windows')
        if audit.get('error') is None:
            raise RuntimeError('No analysis error column available')
        if audit['error'].fillna('').astype(str).str.strip().ne('').any():
            raise RuntimeError('Some candidates failed; raw cache is preserved')
        if audit.get('missing_background_days') is None or (
                audit['missing_background_days'].fillna('')
                .astype(str).str.strip().ne('').any()):
            raise RuntimeError('Missing reference days; raw cache preserved')
        if audit.get('window_context_archive') is None:
            raise RuntimeError('One-second background archives missing')
        for value in audit['window_context_archive']:
            if not isinstance(value, str) or not Path(value).is_file():
                raise RuntimeError('Incomplete one-second context archive')
            validate_gzip(value)
    if len(registry) != len(events):
        raise RuntimeError('Event and registry row counts disagree')
    if events['event_id'].duplicated().any():
        raise RuntimeError('Non-unique event IDs in daily event catalog')
    if set(events['event_id']) != set(registry['event_id']):
        raise RuntimeError('Registry has inconsistent event identities')
    return events, registry, audit


def verify_manifest(manifest):
    """On resume, require the archives and checked outputs still to exist."""
    for key in ('candidate_csv', 'event_csv', 'registry_csv',
                'archive_index', 'context_archive_index'):
        if not Path(manifest[key]).is_file():
            raise RuntimeError('Completed-day file missing: ' + manifest[key])
    for label in ('event', 'registry', 'candidate'):
        path = Path(manifest[f'{label}_csv'])
        if sha256_file(path) != manifest[f'{label}_sha256']:
            raise RuntimeError(f'Checkpoint {label} CSV changed: {path}')
    contexts = pd.read_csv(manifest['context_archive_index'])
    for entry in contexts.to_dict('records'):
        target = Path(entry['archive_file'])
        if not target.is_file() or sha256_file(target) != entry['sha256']:
            raise RuntimeError(f'One-second context archive missing/corrupt: {target}')
        validate_gzip(target)
    index = pd.read_csv(manifest['archive_index'])
    for entry in index.to_dict('records'):
        target = Path(entry['archive_file'])
        if not target.is_file() or sha256_file(target) != entry['sha256']:
            raise RuntimeError(f'Archive verification failed: {target}')
        validate_gzip(target)
    if manifest.get('event_count', -1) != len(index):
        raise RuntimeError('Archived event count differs from checkpoint')
    return True


def summarize(root):
    """Publish a combined overlap-deduplicated catalog after all days finish."""
    files = sorted((root / 'checkpoints').glob('*.json'))
    events_frames, registry_frames, audit_frames, archives = [], [], [], []
    for file in files:
        state = json.loads(file.read_text())
        verify_manifest(state)
        if state['event_count']:
            events_frames.append(pd.read_csv(state['event_csv']))
            registry_frames.append(pd.read_csv(state['registry_csv']))
            archives.append(pd.read_csv(state['archive_index']))
        if state.get('audit_csv') and Path(state['audit_csv']).is_file():
            audit_frame = pd.read_csv(state['audit_csv'])
            audit_frame.insert(0, 'scan_day_utc', state['scan_day_utc'])
            audit_frames.append(audit_frame)
    if events_frames:
        all_events = pd.concat(events_frames, ignore_index=True)
        chosen = pd.DataFrame(deduplicate(all_events.to_dict('records')),
                              columns=FULL_EVENT_COLUMNS)
        selected = set(chosen['event_id'])
        reg = (pd.concat(registry_frames, ignore_index=True)
               .drop_duplicates('event_id'))
        reg = reg.loc[reg.event_id.isin(selected)]
        arc = (pd.concat(archives, ignore_index=True)
               .drop_duplicates('event_id'))
        arc = arc.loc[arc.event_id.isin(selected)]
    else:
        chosen = pd.DataFrame(columns=FULL_EVENT_COLUMNS)
        reg = pd.DataFrame(columns=EVENT_REGISTRY_COLUMNS)
        arc = pd.DataFrame(columns=['event_id','archive_file','sha256'])
    atomic_csv(chosen, root / 'all_events.csv')
    atomic_csv(reg, root / 'all_registry.csv')
    atomic_csv(arc, root / 'all_event_archives.csv')
    atomic_csv(pd.concat(audit_frames, ignore_index=True)
               if audit_frames else pd.DataFrame(),
               root / 'all_audit.csv')
    print(f'[CATALOG] {len(chosen)} distinct candidate events archived',
          flush=True)


def cleanup_raw_cache(start_day, stop_day_exclusive, completed_days,
                      dry_run=False, cache=FULL_CACHE):
    """Remove ONLY canonical local FIELDS full-cadence files on completed days.

    Never recurse, touch user files outside cache, or remove a symlink. The
    caller keeps one preceding processed day during an active scan.
    """
    cache = Path(cache)
    if not cache.is_dir():
        return (0, 0)
    removed = 0
    freed = 0
    for path in sorted(cache.iterdir()):
        if not path.is_file() or path.is_symlink():
            continue
        match = FULL_NAME.fullmatch(path.name)
        if match is None:
            continue
        segment = pd.to_datetime(match.group(1), format='%Y%m%d%H')
        date = segment.normalize()
        if not (start_day <= date < stop_day_exclusive):
            continue
        if date.strftime('%Y%m%d') not in completed_days:
            continue
        if path.parent.resolve() != cache.resolve():
            continue
        length = path.stat().st_size
        if not dry_run:
            path.unlink()
        print(f'[CLEANUP] {"Would remove" if dry_run else "Removed"} '
              f'{path} ({length/1048576:.1f} MiB)', flush=True)
        removed += 1
        freed += length
    return removed, freed


def run_day(root, start, end, day, position_source='cached-spc',
            context_minutes=5, plots=False):
    key = day.strftime('%Y%m%d')
    name = f'{root.name}_{key}_candidates.csv'
    catalog = root / 'candidates' / name
    catalog.parent.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, 'code/magnetic_data_reader.py',
           '--scan', '--start', start.isoformat(), '--end', end.isoformat(),
           '--output', str(catalog)]
    print(f'[PIPELINE] Reader {start} to {end} UTC', flush=True)
    subprocess.run(cmd, check=True)

    # The reader writes actual one-second coverage to a separate checkpoint.
    # A sparse day can be analyzed but MUST NOT authorize raw-cache deletion.
    scan_checkpoint = Path(str(catalog) + '.checkpoint.json')
    observed_fraction = 0.
    if scan_checkpoint.is_file():
        scan_status = json.loads(scan_checkpoint.read_text())
        expected = max(1, int(np.ceil((end - start).total_seconds())))
        observed_fraction = int(scan_status.get('valid_seconds', 0)) / expected
    cleanup_eligible = observed_fraction >= .90
    if not cleanup_eligible:
        print(f'[COVERAGE] {day.date()}: only '
              f'{observed_fraction:.1%} valid 1-second bins. '
              'Will not delete raw files from this day.', flush=True)

    events_path, registry_path, audit_path = analysis_outputs(catalog)
    seeds = pd.read_csv(catalog)
    if not seeds.empty:
        cmd = [sys.executable, 'code/magnetic_data_analysis.py',
               '--windows-csv', str(catalog), '--position-source',
               position_source, '--window-archive-dir',
               str(root / 'contexts' / key)]
        if not plots:
            cmd.append('--no-plots')
        print(f'[PIPELINE] Analysis: {len(seeds)} candidate windows',
              flush=True)
        subprocess.run(cmd, check=True)
    else:
        # A valid quiet day must still have complete empty result artifacts.
        atomic_csv(pd.DataFrame(columns=FULL_EVENT_COLUMNS), events_path)
        atomic_csv(pd.DataFrame(columns=EVENT_REGISTRY_COLUMNS), registry_path)
        atomic_csv(pd.DataFrame(columns=['seed_utc', 'error']), audit_path)
    events, registry, audit = validate_analysis(
        catalog, events_path, registry_path, audit_path)

    archive_dir = root / 'native_events' / key
    archives = []
    for i, row in enumerate(events.to_dict('records'), 1):
        print(f'[ARCHIVE {i}/{len(events)}] {row["event_id"]}',
              flush=True)
        archives.append(archive_native_event(row, archive_dir, context_minutes))
    index_path = root / 'archive_indices' / f'{key}.csv'
    context_index_path = root / 'context_indices' / f'{key}.csv'
    contexts = []
    if not audit.empty:
        for source in audit['window_context_archive'].tolist():
            contexts.append(dict(archive_file=str(source),
                                 sha256=sha256_file(source)))
    atomic_csv(pd.DataFrame(contexts, columns=['archive_file','sha256']),
               context_index_path)
    atomic_csv(pd.DataFrame(archives, columns=[
        'event_id','archive_file','first_requested_utc',
        'end_requested_exclusive_utc','native_archived_samples',
        'native_event_archived_samples','native_event_valid_samples','sha256'
    ]), index_path)
    state = dict(scan_day_utc=day.strftime('%Y-%m-%d'),
                 start_utc=start.isoformat(), end_utc=end.isoformat(),
                 candidate_csv=str(catalog), event_csv=str(events_path),
                 registry_csv=str(registry_path),
                 candidate_sha256=sha256_file(catalog),
                 event_sha256=sha256_file(events_path),
                 registry_sha256=sha256_file(registry_path),
                 audit_csv=str(audit_path) if not seeds.empty else '',
                 context_archive_index=str(context_index_path),
                 archive_index=str(index_path), event_count=len(events),
                 screening_windows=len(seeds),
                 valid_seconds_coverage_fraction=round(observed_fraction, 6),
                 cleanup_eligible=cleanup_eligible,
                 native_sample_archive='quality_filtered_RTNs_at_native_cadence',
                 context_minutes=context_minutes)
    verify_manifest(state)
    checkpoint = root / 'checkpoints' / f'{key}.json'
    atomic_json(state, checkpoint)
    return state


def pipeline(args):
    begin = _as_utc_naive(args.start)
    stop = (_parse_duration(begin, args.duration)
            if args.duration else _as_utc_naive(args.end))
    if stop <= begin:
        raise ValueError('Invalid requested observation interval')
    root = DATA / f'rolling_{begin:%Y%m%dT%H%M}_{stop:%Y%m%dT%H%M}'
    first = begin.normalize()
    last = (stop - pd.Timedelta(nanoseconds=1)).normalize()
    days = pd.date_range(first, last, freq='D')
    if root.exists() and not args.resume and (
            list((root / 'checkpoints').glob('*.json')) if
            (root / 'checkpoints').exists() else False):
        raise RuntimeError('This pipeline already has saved progress. '
                           'Use --resume to protect existing archives.')
    root.mkdir(parents=True, exist_ok=True)
    print(f'[PIPELINE] {begin} to {stop}, {len(days)} UTC days, '
          f'cleanup={"ON" if args.cleanup_raw else "OFF"}', flush=True)
    completed = set()
    cleanup_eligible_days = set()
    for i, day in enumerate(days, 1):
        key = day.strftime('%Y%m%d')
        checkpoint = root / 'checkpoints' / f'{key}.json'
        if checkpoint.exists():
            if not args.resume:
                raise RuntimeError('Existing completed day; use --resume')
            state = json.loads(checkpoint.read_text())
            verify_manifest(state)
            print(f'[RESUME] {key}: verified archived day; skipping',
                  flush=True)
        else:
            start = max(begin, day)
            end = min(stop, day + pd.Timedelta(days=1))
            print(f'[DAY {i}/{len(days)}] {key}', flush=True)
            state = run_day(root, start, end, day,
                            position_source=args.position_source,
                            context_minutes=args.context_minutes,
                            plots=args.plots)
        completed.add(key)
        if state.get('cleanup_eligible', False):
            cleanup_eligible_days.add(key)
        if args.cleanup_raw:
            # Keep the previous day available for the next day's 100-minute
            # local baseline and for event windows spanning midnight.
            safe_stop = day - pd.Timedelta(days=1)
            count, total = cleanup_raw_cache(
                first, safe_stop, cleanup_eligible_days,
                dry_run=args.dry_run_cleanup)
            if count:
                print(f'[CLEANUP] {count} CDFs; '
                      f'{total/1073741824:.2f} GiB '
                      f'{"eligible" if args.dry_run_cleanup else "reclaimed"}',
                      flush=True)
    # A combined catalog is built only once. Daily catalogs/checkpoints
    # already exist after every day; this avoids quadratic repeated
    # deduplication of thousands of events during a long scan.
    summarize(root)
    if args.cleanup_raw:
        count, total = cleanup_raw_cache(
            first, last + pd.Timedelta(days=1), cleanup_eligible_days,
            dry_run=args.dry_run_cleanup)
        print(f'[FINAL CLEANUP] {count} additional CDFs; '
              f'{total/1073741824:.2f} GiB '
              f'{"eligible" if args.dry_run_cleanup else "reclaimed"}',
              flush=True)
    print(f'[COMPLETE] Consolidated catalogs and verified archives: {root}',
          flush=True)


def self_test():
    with TemporaryDirectory() as directory:
        home = Path(directory)
        cache = home / 'fields_full'
        cache.mkdir()
        for date in ('2018110100', '2018110200', '2018110300'):
            (cache / f'psp_fld_l2_mag_rtn_{date}_v02.cdf').write_bytes(
                b'Synthetic CDF placeholder: never read')
        (cache / 'unrelated.txt').write_text('Must keep')
        completed = {'20181101', '20181102'}
        first = pd.Timestamp('2018-11-01')
        count, _ = cleanup_raw_cache(
            first, pd.Timestamp('2018-11-02'), completed,
            cache=cache, dry_run=True)
        assert count == 1
        assert len(list(cache.glob('*.cdf'))) == 3
        count, _ = cleanup_raw_cache(
            first, pd.Timestamp('2018-11-03'), completed, cache=cache)
        assert count == 2
        assert len(list(cache.glob('*.cdf'))) == 1
        assert (cache / 'unrelated.txt').is_file()

        idx = pd.date_range('2018-11-02 12:00:00', periods=20, freq='s')
        raw = pd.DataFrame({'BR':1., 'BT':0., 'BN':0., 'Bmag':1.},
                           index=idx)
        event = dict(event_id='PSP-MAG-TEST',
                     event_start_utc=idx[5].isoformat(),
                     event_end_utc=idx[14].isoformat(),
                     native_valid_samples=10)
        archive = archive_native_event(
            event, home / 'archived', context_minutes=0,
            stream_reader=lambda left,right: iter([raw]))
        with gzip.open(archive['archive_file'], 'rt') as handle:
            saved = pd.read_csv(handle)
        assert len(saved) == 10
        assert saved.iloc[0]['UTC'] == idx[5].isoformat(sep=' ')
        assert sha256_file(archive['archive_file']) == archive['sha256']
    print('PASS: archive sample integrity, cache scoping and safe cleanup')


def main():
    parser = argparse.ArgumentParser(
        description='Daily FIELDS reader -> analysis -> native event archives '
                    '-> optional verified cache cleanup')
    parser.add_argument('--start', default='2018-11-01',
                        help='UTC start date or datetime')
    period = parser.add_mutually_exclusive_group()
    period.add_argument('--duration', default=None,
                        help='e.g. 2mo, 10d, 3w')
    period.add_argument('--end', default=None,
                        help='UTC exclusive end')
    parser.add_argument('--cleanup-raw', action='store_true',
                        help='Delete old 6-hour CDF cache files only after '
                             'all candidate samples have been archived')
    parser.add_argument('--dry-run-cleanup', action='store_true',
                        help='Report eligible CDF files without deleting them')
    parser.add_argument('--resume', action='store_true',
                        help='Verify existing completed days and continue')
    parser.add_argument('--position-source',
                        choices=('cached-spc','spc','none'),
                        default='cached-spc')
    parser.add_argument('--context-minutes', type=int, default=5,
                        help='Native original-cadence samples to preserve '
                             'before and after each event')
    parser.add_argument('--plots', action='store_true',
                        help='Write analysis PNGs (off by default)')
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if not args.duration and not args.end:
        parser.error('Specify --duration or --end')
    if args.context_minutes < 0:
        parser.error('--context-minutes must be nonnegative')
    if args.dry_run_cleanup and not args.cleanup_raw:
        parser.error('--dry-run-cleanup requires --cleanup-raw')
    pipeline(args)


if __name__ == '__main__':
    main()
