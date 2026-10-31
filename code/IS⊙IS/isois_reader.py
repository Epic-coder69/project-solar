#!/usr/bin/env python3
"""Parker Solar Probe IS⊙IS unified reader.

Contains the *actual* L2 SUMMARY and EPI-Hi CDF decoding implementations.
Preserves source CDF integrations, energy bins, metadata and missing values.
This module performs ingestion only; analysis belongs to isois_analyzer.py.
Previous cache and output paths remain unchanged.
"""
from __future__ import annotations

# === L2 SUMMARY CDF READER (previously particle_data_reader.py) ===
import argparse
import hashlib
import json
import re
import sys
import time
from datetime import date, timedelta
from pathlib import Path
from urllib.parse import urljoin

import numpy as np
import pandas as pd

BASE_URL = 'https://spp-isois.sr.unh.edu/data_public/ISOIS/level2/'
PRODUCT = 'PSP_ISOIS_L2-SUMMARY'
FILE_PATTERN = re.compile(r'^psp_isois_l2-summary_(\d{8})_v(\d+)\.cdf$')
DEFAULT_VARIABLES = (
    'A_H_Rate_TS', 'A_Heavy_Rate_TS',
    'HET_A_H_Rate_TS', 'HET_A_Electrons_Rate_TS',
)
DEFAULT_CACHE = Path('data/isois/summary_cdf')
OUTPUT_COLUMNS = (
    'timestamp_utc', 'integration_start_utc', 'integration_end_utc',
    'variable', 'value', 'units', 'valid_measurement', 'quality_reason',
    'source_file', 'source_version', 'product', 'analysis_grade',
)
EVENT_WINDOW_COLUMNS = (
    'event_id', 'event_start_utc', 'event_end_exclusive_utc',
    'window_start_utc', 'window_end_exclusive_utc', 'boundary_basis',
)
SCHEMA_VERSION = 2


def _cdf_module():
    try:
        import cdflib
    except ImportError as exc:
        raise RuntimeError('Missing cdflib: python -m pip install -r code/FIELDS/external_libraries_requirements.txt') from exc
    return cdflib


def _sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def _atomic_text(path, content):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.tmp')
    try:
        temp.write_text(content, encoding='utf-8')
        temp.replace(path)
    finally:
        temp.unlink(missing_ok=True)


def _atomic_csv(frame, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.tmp')
    try:
        frame.to_csv(temp, index=False)
        temp.replace(path)
    finally:
        temp.unlink(missing_ok=True)


def _attr(attrs, name, default=None):
    value = attrs.get(name, default)
    if isinstance(value, np.ndarray) and value.size == 1:
        value = value.reshape(-1)[0]
    if isinstance(value, (tuple, list)) and len(value) == 1:
        value = value[0]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, bytes):
        value = value.decode('utf-8', errors='replace')
    return value


def _number(value):
    if value is None:
        return None
    try:
        value = float(np.asarray(value).reshape(-1)[0])
        return value if np.isfinite(value) else None
    except (ValueError, TypeError, IndexError):
        return None


def _utc(value):
    result = pd.to_datetime(value, utc=True, errors='raise')
    if pd.isna(result):
        raise ValueError('Missing UTC timestamp')
    return result


def _iso(timestamp):
    return timestamp.isoformat().replace('+00:00', 'Z')


def _epoch_nanoseconds(series):
    """Normalize UTC timestamps to epoch ns across pandas 2.x/3.x resolutions."""
    return pd.to_datetime(series, utc=True, errors='raise').dt.as_unit('ns').astype('int64')


def _half_width_seconds(raw, attrs):
    """Fail closed on ambiguous time units (do not mistake ns for s)."""
    units = str(_attr(attrs, 'UNITS', '')).lower().strip()
    si = str(_attr(attrs, 'SI_conversion', '')).lower().replace(' ', '')
    if units in ('ns', 'nsec', 'nanosecond', 'nanoseconds') or si.startswith('1e-9>'):
        multiplier = 1e-9
    elif units in ('ms', 'msec', 'millisecond', 'milliseconds') or si.startswith('1e-3>'):
        multiplier = 1e-3
    elif units in ('s', 'sec', 'second', 'seconds') or si.startswith('1>'):
        multiplier = 1.0
    else:
        raise ValueError(f'Unrecognized integration delta unit {units!r}, SI_conversion {si!r}')
    seconds = np.asarray(raw, dtype=float) * multiplier
    if seconds.ndim != 1 or not np.isfinite(seconds).all() or (seconds <= 0).any():
        raise ValueError('Nonpositive, nonfinite, or multidimensional integration half widths')
    return seconds


def extract_rate(cdf, variable, source_filename, *, epoch_converter):
    """Decode one 1D survey count-rate variable, including integration bounds."""
    attributes = cdf.varattsget(variable)
    epoch_var = _attr(attributes, 'DEPEND_0')
    units = str(_attr(attributes, 'UNITS', '')).strip()
    if not epoch_var or not units:
        raise ValueError(f'{variable}: missing DEPEND_0 or UNITS')
    if 'rate' not in variable.lower() or ('count' not in units.lower() and 's!' not in units.lower()):
        raise ValueError(f'{variable}: unsupported rate metadata ({units})')
    values = np.asarray(cdf.varget(variable))
    if values.ndim != 1:
        raise ValueError(f'{variable}: expected 1D rates; found shape {values.shape}; refuse to collapse bins')
    epoch_raw = cdf.varget(str(epoch_var))
    epoch = pd.DatetimeIndex(pd.to_datetime(epoch_converter(epoch_raw), utc=True, errors='coerce'))
    if len(epoch) != len(values) or epoch.isna().any() or epoch.has_duplicates or not epoch.is_monotonic_increasing:
        raise ValueError(f'{variable}: missing, nonmonotonic, duplicate, or unmatched epoch tags')
    epoch_attrs = cdf.varattsget(str(epoch_var))
    delta_names = [(_attr(epoch_attrs, key), key) for key in ('DELTA_MINUS_VAR', 'DELTA_PLUS_VAR')]
    if any(not name for name, _ in delta_names):
        raise ValueError(f'{variable}: missing DELTA_PLUS_VAR/DELTA_MINUS_VAR on {epoch_var}')
    halfwidths = []
    for name, _ in delta_names:
        raw = np.asarray(cdf.varget(str(name)))
        if raw.ndim == 0:
            raw = np.repeat(raw, len(values))
        if raw.shape != values.shape:
            raise ValueError(f'{name}: delta must be scalar or one value per sample')
        halfwidths.append(_half_width_seconds(raw, cdf.varattsget(str(name))))
    start = epoch - pd.to_timedelta(halfwidths[0], unit='s')
    end = epoch + pd.to_timedelta(halfwidths[1], unit='s')
    if start.isna().any() or end.isna().any():
        raise ValueError(f'{variable}: invalid integration timestamps')
    data = pd.to_numeric(pd.Series(values), errors='coerce').to_numpy(dtype=float, copy=True)
    invalid = ~np.isfinite(data)
    reasons = np.full(len(data), '', dtype=object)
    reasons[invalid] = 'nonfinite'
    fillval = _number(_attr(attributes, 'FILLVAL'))
    if fillval is not None:
        mask = np.isclose(data, fillval, rtol=1e-7, atol=0) & ~invalid
        invalid |= mask
        reasons[mask] = 'fill_value'
    minimum = _number(_attr(attributes, 'VALIDMIN'))
    maximum = _number(_attr(attributes, 'VALIDMAX'))
    for threshold, cause, condition in (
        (minimum, 'below_valid_min', lambda x, y: x < y),
        (maximum, 'above_valid_max', lambda x, y: x > y),
        (0, 'negative_rate', lambda x, y: x < y),
    ):
        if threshold is not None:
            mask = condition(data, threshold) & np.isfinite(data) & ~invalid
            invalid |= mask
            reasons[mask] = cause
    data[invalid] = np.nan
    match = FILE_PATTERN.fullmatch(Path(source_filename).name)
    if not match:
        raise ValueError(f'Unexpected summary CDF filename: {source_filename}')
    return pd.DataFrame({
        'timestamp_utc': [_iso(x) for x in epoch],
        'integration_start_utc': [_iso(x) for x in start],
        'integration_end_utc': [_iso(x) for x in end],
        'variable': variable,
        'value': data,
        'units': units,
        'valid_measurement': ~invalid,
        'quality_reason': reasons,
        'source_file': Path(source_filename).name,
        'source_version': int(match.group(2)),
        'product': PRODUCT,
        'analysis_grade': 'survey_only',
    }, columns=OUTPUT_COLUMNS)


def read_cdf(path, variables=DEFAULT_VARIABLES):
    """Open CDF once; reject missing channels, bad dimensions and timestamps."""
    cdf_module = _cdf_module()
    path = Path(path)
    if not FILE_PATTERN.fullmatch(path.name):
        raise ValueError(f'Not an ISʘIS summary CDF name: {path.name}')
    cdf = cdf_module.CDF(str(path))
    try:
        info = cdf.cdf_info()
        available = set(info.zVariables) | set(info.rVariables)
        missing = set(variables) - available
        if missing:
            raise ValueError(f'{path.name}: missing requested channels: {sorted(missing)}')
        frames = [extract_rate(cdf, var, path.name, epoch_converter=cdf_module.cdfepoch.to_datetime)
                  for var in variables]
    finally:
        if hasattr(cdf, 'close'):
            cdf.close()
    if not frames:
        raise ValueError('No channels requested')
    return pd.concat(frames, ignore_index=True)


def inspect_file(path):
    cdf_module = _cdf_module()
    cdf = cdf_module.CDF(str(path))
    try:
        info = cdf.cdf_info()
        for name in list(info.zVariables) + list(info.rVariables):
            attrs = cdf.varattsget(name)
            print(f'{name:44s} epoch={str(_attr(attrs, "DEPEND_0", "")):18s} '
                  f'units={str(_attr(attrs, "UNITS", "")):23s} '
                  f'{str(_attr(attrs, "CATDESC", ""))[:55]}')
    finally:
        if hasattr(cdf, 'close'):
            cdf.close()


def validate_fields_events(path, *, allow_unverified=False):
    """Read FIELDS strong events or registry, checking its retention manifest.

    FIELDS raw event end is inclusive at 1 s. Registry end_exclusive is
    authoritative; convert raw event end by +1 s, never silently by +0 s.
    """
    path = Path(path)
    if not (path.name.endswith('_strong_events.csv') or path.name.endswith('_strong_registry.csv')):
        raise ValueError('Use a FIELDS *_strong_events.csv or *_strong_registry.csv')
    if not path.is_file():
        raise FileNotFoundError(path)
    manifest_path = path.with_name(path.name.replace('_strong_events.csv', '_strong_retention.json')
                                   .replace('_strong_registry.csv', '_strong_retention.json'))
    if not allow_unverified:
        if not manifest_path.is_file():
            raise ValueError(f'Missing FIELDS retention manifest {manifest_path}')
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        check_key = ('strong_registry_sha256' if path.name.endswith('_strong_registry.csv')
                     else 'strong_events_sha256')
        if manifest.get(check_key) != _sha256(path):
            raise ValueError(f'FIELDS catalog integrity check failed: {check_key}')
    else:
        manifest = {}
    frame = pd.read_csv(path)
    required = {'event_id', 'event_start_utc', 'event_end_utc'}
    missing = required - set(frame.columns)
    if missing or frame.empty:
        raise ValueError(f'Empty or invalid FIELDS event catalog; missing={sorted(missing)}')
    if frame['event_id'].isna().any() or frame['event_id'].duplicated().any():
        raise ValueError('Missing or duplicated FIELDS event IDs')
    if manifest.get('strong_events_retained') is not None and len(frame) != int(manifest['strong_events_retained']):
        raise ValueError('FIELDS event count differs from retention manifest')
    rows = []
    registry = path.name.endswith('_strong_registry.csv')
    for item in frame.to_dict('records'):
        start = _utc(item['event_start_utc'])
        if registry:
            if 'event_end_exclusive_utc' not in item or pd.isna(item['event_end_exclusive_utc']):
                raise ValueError('Registry missing its authoritative event_end_exclusive_utc')
            end = _utc(item['event_end_exclusive_utc'])
            if end != _utc(item['event_end_utc']) + pd.Timedelta(seconds=1):
                raise ValueError('FIELDS registry inclusive/exclusive ends are inconsistent')
            basis = 'fields_registry_exclusive'
        else:
            end = _utc(item['event_end_utc']) + pd.Timedelta(seconds=1)
            basis = 'fields_events_inclusive_plus_one_second'
        if end <= start:
            raise ValueError(f'Nonpositive FIELDS event interval: {item["event_id"]}')
        rows.append((str(item['event_id']), start, end, basis))
    return rows


def build_event_windows(events, *, before_minutes=30, after_minutes=30):
    if not (np.isfinite(before_minutes) and np.isfinite(after_minutes)) or min(before_minutes, after_minutes) < 0:
        raise ValueError('Padding must be nonnegative and finite')
    rows = []
    for event_id, start, end, basis in events:
        first = start - pd.Timedelta(minutes=before_minutes)
        last = end + pd.Timedelta(minutes=after_minutes)
        rows.append({'event_id': event_id, 'event_start_utc': _iso(start),
                     'event_end_exclusive_utc': _iso(end), 'window_start_utc': _iso(first),
                     'window_end_exclusive_utc': _iso(last), 'boundary_basis': basis})
    return pd.DataFrame(rows, columns=EVENT_WINDOW_COLUMNS)


def merged_intervals(windows):
    """Nonoverlapping half-open windows; adjacent windows merge too."""
    intervals = sorted((_utc(r.window_start_utc), _utc(r.window_end_exclusive_utc))
                       for r in windows.itertuples(index=False))
    merged = []
    for start, end in intervals:
        if end <= start:
            raise ValueError('Invalid event observation interval')
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


def required_days(intervals, *, halo_days=1):
    """Daily ISʘIS files are organized by midpoint timestamp.

    Include neighbor days so an integration whose midpoint lies across a
    UTC boundary is not lost. The conservatism costs small summary files.
    """
    if halo_days not in (0, 1):
        raise ValueError('Day halo must be 0 or 1')
    result = set()
    for start, end in intervals:
        a = start.date() - timedelta(days=halo_days)
        last = (end - pd.Timedelta(nanoseconds=1)).date() + timedelta(days=halo_days)
        while a <= last:
            result.add(a)
            a += timedelta(days=1)
    return sorted(result)


def filter_by_overlap(frame, intervals):
    """Return unique integrations that intersect any requested [start,end).

    Uses time intervals, NOT the midpoint. A single integration may overlap
    several FIELDS events but is stored only once in the shared rate table.
    """
    if frame.empty or not intervals:
        return frame.iloc[:0].copy()
    left = _epoch_nanoseconds(frame['integration_start_utc']).to_numpy()
    right = _epoch_nanoseconds(frame['integration_end_utc']).to_numpy()
    if (right <= left).any():
        raise ValueError('Bad integration interval in decoded CDF')
    starts = np.asarray([a.value for a, _ in intervals], dtype=np.int64)
    ends = np.asarray([b.value for _, b in intervals], dtype=np.int64)
    index = np.searchsorted(ends, left, side='right')
    inside = index < len(intervals)
    eligible = np.flatnonzero(inside)
    inside[eligible] = starts[index[eligible]] < right[eligible]
    return frame.iloc[np.flatnonzero(inside)].copy()


def source_index(html):
    versions = {}
    for day, ver in re.findall(r'psp_isois_l2-summary_(\d{8})_v(\d+)\.cdf', html, flags=re.I):
        name = f'psp_isois_l2-summary_{day}_v{ver}.cdf'
        if day not in versions or int(ver) > versions[day][0]:
            versions[day] = (int(ver), name)
    return {k: v[1] for k, v in versions.items()}


def choose_sources(days, cache, *, download, offline, session=None):
    """Prefer highest version advertised; permit explicitly offline cache."""
    cache = Path(cache)
    if offline and download:
        raise ValueError('--offline and --download cannot both be set')
    remote = {}
    if not offline:
        if session is None:
            raise ValueError('A requests session is required for online discovery')
        r = session.get(BASE_URL, timeout=(15, 90))
        r.raise_for_status()
        remote = source_index(r.text)
        if not remote:
            raise RuntimeError('No summary-product filenames found in the server index')
    local = {}
    for p in cache.glob('psp_isois_l2-summary_*_v*.cdf'):
        match = FILE_PATTERN.fullmatch(p.name)
        if match and p.is_file() and p.stat().st_size:
            key, version = match.group(1), int(match.group(2))
            if key not in local or version > local[key][0]:
                local[key] = (version, p.name)
    selected, missing = {}, []
    for day in days:
        key = day.strftime('%Y%m%d')
        chosen = (local[key][1] if offline and key in local
                  else remote.get(key) if not offline else None)
        if chosen is None:
            missing.append(day.isoformat())
        elif not (cache / chosen).exists() and not download:
            missing.append(day.isoformat())
        else:
            selected[day] = chosen
    return selected, missing


def download_file(session, filename, cache, *, attempts=4, pause=time.sleep):
    """Verified atomic HTTP transfer with safe Range resume and no deletion."""
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    dest = cache / filename
    part = cache / (filename + '.part')
    if dest.is_file() and dest.stat().st_size > 0:
        return dest
    for trial in range(1, attempts + 1):
        offset = part.stat().st_size if part.exists() else 0
        headers = {'Range': f'bytes={offset}-'} if offset else {}
        try:
            with session.get(urljoin(BASE_URL, filename), stream=True,
                             timeout=(20, 120), headers=headers) as response:
                response.raise_for_status()
                if offset and response.status_code == 206:
                    match = re.match(r'^bytes (\d+)-(\d+)/(\d+)$',
                                     response.headers.get('Content-Range', ''))
                    if not match or int(match.group(1)) != offset:
                        raise IOError('Server sent invalid Content-Range; cannot safely append')
                    mode, expected_total = 'ab', int(match.group(3))
                else:
                    # Server ignored Range (HTTP 200) or this is a fresh transfer.
                    if response.status_code != 200:
                        raise IOError(f'Unexpected status {response.status_code}')
                    mode = 'wb'
                    expected_total = int(response.headers.get('Content-Length') or 0)
                with part.open(mode) as sink:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            sink.write(chunk)
                size = part.stat().st_size
                if size == 0 or expected_total and size != expected_total:
                    raise IOError(f'Transfer incomplete: {size} bytes, expected {expected_total}')
            # Open the CDF before making it the trusted cached file.
            module = _cdf_module()
            opened = module.CDF(str(part))
            try:
                info = opened.cdf_info()
                if not info.zVariables and not info.rVariables:
                    raise IOError('Downloaded CDF has no variables')
            finally:
                if hasattr(opened, 'close'):
                    opened.close()
            part.replace(dest)
            return dest
        except Exception as exc:
            print(f'[DOWNLOAD RETRY {trial}/{attempts}] {filename}: {exc}', file=sys.stderr)
            if trial >= attempts:
                raise
            pause(min(2 * trial, 8))
    raise AssertionError('Unreachable')


def _partition_paths(output, day):
    folder = output.parent / (output.stem + '.daily_parts')
    return folder / (day.strftime('%Y%m%d') + '.csv'), folder / (day.strftime('%Y%m%d') + '.json')


def _fingerprint(*, variables, intervals, catalog_hash, day_halo):
    payload = {'variables': list(variables), 'intervals': [[_iso(a), _iso(b)] for a, b in intervals],
               'fields_catalog_sha256': catalog_hash, 'day_halo': day_halo,
               'reader_schema_version': SCHEMA_VERSION}
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode('utf-8')).hexdigest()


def write_partition(frame, output, day, source_file, fingerprint):
    csv, metadata = _partition_paths(output, day)
    _atomic_csv(frame, csv)
    record = {'day': day.isoformat(), 'cdf_file': source_file.name,
              'cdf_sha256': _sha256(source_file), 'csv_sha256': _sha256(csv),
              'rows': len(frame), 'fingerprint': fingerprint}
    _atomic_text(metadata, json.dumps(record, indent=2) + '\n')
    return csv


def reusable_partition(output, day, source_file, fingerprint):
    csv, meta = _partition_paths(output, day)
    if not csv.is_file() or not meta.is_file():
        return None
    try:
        record = json.loads(meta.read_text(encoding='utf-8'))
        if (record.get('fingerprint') == fingerprint
                and record.get('cdf_file') == source_file.name
                and record.get('cdf_sha256') == _sha256(source_file)
                and record.get('csv_sha256') == _sha256(csv)):
            if not set(OUTPUT_COLUMNS).issubset(pd.read_csv(csv, nrows=0).columns):
                return None
            return csv
    except (OSError, ValueError, KeyError):
        return None
    return None



def coverage_audit(windows, measurements, variables):
    """Coverage only, NOT particle-event detection or a statistical test.

    Union overlapping integration intervals to avoid overcounting exposure.
    This gives explicit pre/during/post availability for every FIELDS event.
    """
    by_channel = {}
    for variable, group in measurements.groupby('variable'):
        group = group.copy()
        group['_start_ns'] = _epoch_nanoseconds(group.integration_start_utc)
        group['_end_ns'] = _epoch_nanoseconds(group.integration_end_utc)
        group = group.sort_values('_end_ns')
        starts = group['_start_ns'].to_numpy(dtype=np.int64)
        ends = group['_end_ns'].to_numpy(dtype=np.int64)
        values = pd.to_numeric(group.value, errors='coerce').to_numpy(dtype=float)
        valid = group.valid_measurement.astype(str).str.lower().eq('true').to_numpy() & np.isfinite(values)
        if (ends <= starts).any():
            raise ValueError(f'{variable}: invalid source integration intervals')
        by_channel[variable] = (starts, ends, valid, int(np.max(ends - starts)))
    audit = []
    for w in windows.itertuples(index=False):
        t0 = _utc(w.window_start_utc)
        a = _utc(w.event_start_utc)
        b = _utc(w.event_end_exclusive_utc)
        t1 = _utc(w.window_end_exclusive_utc)
        for label, left, right in (('before', t0, a), ('during', a, b), ('after', b, t1)):
            if left == right:
                continue
            left_ns, right_ns = left.value, right.value
            duration = (right_ns - left_ns) * 1e-9
            for variable in variables:
                if variable in by_channel:
                    starts, ends, good, max_width = by_channel[variable]
                    lo = int(np.searchsorted(ends, left_ns, side='right'))
                    hi = int(np.searchsorted(ends, right_ns + max_width, side='left'))
                    candidates = (starts[lo:hi] < right_ns) & (ends[lo:hi] > left_ns)
                    count = int(np.count_nonzero(candidates))
                    valid_indices = np.flatnonzero(candidates & good[lo:hi]) + lo
                    spans = sorted((max(int(starts[i]), left_ns), min(int(ends[i]), right_ns))
                                   for i in valid_indices)
                    union_ns = 0
                    if spans:
                        begin, finish = spans[0]
                        for low, high in spans[1:]:
                            if low <= finish:
                                finish = max(finish, high)
                            else:
                                union_ns += finish - begin
                                begin, finish = low, high
                        union_ns += finish - begin
                    exposure = union_ns * 1e-9
                else:
                    count, valid_indices, exposure = 0, [], 0.0
                coverage = min(1.0, exposure / duration)
                audit.append({
                    'event_id': w.event_id, 'variable': variable, 'segment': label,
                    'segment_start_utc': _iso(left), 'segment_end_exclusive_utc': _iso(right),
                    'valid_integration_seconds': round(exposure, 6),
                    'coverage_fraction': round(coverage, 6),
                    'overlapping_samples': count,
                    'valid_samples': int(len(valid_indices)),
                    'status': ('no_measurement' if not count else
                               'no_valid_samples' if not len(valid_indices) else
                               'partial_coverage' if coverage < .999999 else 'covered'),
                })
    return pd.DataFrame(audit)


def event_mode(args):
    events = validate_fields_events(args.events_csv, allow_unverified=args.allow_unverified_events)
    windows = build_event_windows(events, before_minutes=args.before_minutes,
                                  after_minutes=args.after_minutes)
    intervals = merged_intervals(windows)
    days = required_days(intervals, halo_days=args.day_halo)
    output = args.output or Path('data/isois/fields_event_summary_rates.csv')
    print(f'[PLAN] {len(windows)} FIELDS events -> {len(intervals)} merged observation windows; '
          f'{len(days)} source UTC days with a {args.day_halo}-day boundary halo')
    if args.plan_only:
        for left, right in intervals[:8]:
            print(f'  {_iso(left)} -> {_iso(right)}')
        print('[PLAN ONLY] No NASA access, data download, or output changes')
        return
    window_path = output.with_suffix('.event_windows.csv')
    _atomic_csv(windows, window_path)
    catalog_hash = _sha256(args.events_csv)
    fingerprint = _fingerprint(variables=args.variables, intervals=intervals,
                               catalog_hash=catalog_hash, day_halo=args.day_halo)
    process_sources(args, days=days, intervals=intervals, output=output,
                    fingerprint=fingerprint, catalog_hash=catalog_hash,
                    event_count=len(events), merged_window_count=len(intervals),
                    event_windows=windows)


def process_sources(args, *, days, intervals, output, fingerprint, catalog_hash=None,
                    event_count=None, merged_window_count=None, event_windows=None):
    import requests
    with requests.Session() as session:
        selected, missing = choose_sources(days, args.cache, download=args.download,
                                           offline=args.offline, session=session)
        partition_files, source_records = [], []
        for index, (day, filename) in enumerate(selected.items(), 1):
            path = args.cache / filename
            if not path.exists():
                path = download_file(session, filename, args.cache)
            if not path.is_file() or not path.stat().st_size:
                raise RuntimeError(f'Empty or missing CDF cache {path}')
            if args.resume:
                reused = reusable_partition(output, day, path, fingerprint)
                if reused is not None:
                    print(f'[RESUME {index}/{len(selected)}] {day}: verified partition reused')
                    part = reused
                else:
                    part = None
            else:
                part = None
            if part is None:
                decoded = read_cdf(path, args.variables)
                if intervals is not None:
                    decoded = filter_by_overlap(decoded, intervals)
                else:
                    from_day = pd.Timestamp(days[0], tz='UTC')
                    until_day = pd.Timestamp(days[-1] + timedelta(days=1), tz='UTC')
                    mid = pd.to_datetime(decoded.timestamp_utc, utc=True)
                    decoded = decoded.loc[(mid >= from_day) & (mid < until_day)].copy()
                part = write_partition(decoded, output, day, path, fingerprint)
                print(f'[READ {index}/{len(selected)}] {day}: {len(decoded)} selected rows')
            partition_files.append(part)
            record = json.loads(_partition_paths(output, day)[1].read_text(encoding='utf-8'))
            source_records.append({'utc_file_day': day.isoformat(), 'filename': filename,
                                   'source_sha256': record['cdf_sha256'], 'selected_rows': int(record['rows'])})
    if not partition_files:
        raise RuntimeError('No source days available; no output was overwritten')
    tables = [pd.read_csv(p) for p in partition_files]
    nonempty = [table for table in tables if not table.empty]
    if not nonempty:
        raise RuntimeError('Zero selected measurements: no output written')
    final = pd.concat(nonempty, ignore_index=True)
    if final.empty:
        raise RuntimeError('Zero selected measurements: no output written')
    mid = pd.to_datetime(final['timestamp_utc'], utc=True)
    final['_sort'] = mid
    final = final.sort_values(['variable', '_sort']).drop(columns=['_sort']).reset_index(drop=True)
    duplicate = final.duplicated(['variable', 'timestamp_utc'])
    if duplicate.any():
        raise ValueError(f'{int(duplicate.sum())} duplicate channel/epoch records; no output overwritten')
    if not set(final['product']) == {PRODUCT} or not set(final['analysis_grade']) == {'survey_only'}:
        raise ValueError('Unexpected mixed particle product; no output overwritten')
    # Audit per-event coverage without changing or reclassifying FIELDS events.
    audit_file = None
    if event_windows is not None:
        audit = coverage_audit(event_windows, final, args.variables)
        audit_file = output.with_suffix('.coverage_audit.csv')
        _atomic_csv(audit, audit_file)
        print(f'[COVERAGE] {len(audit)} event/channel/segment availability rows: {audit_file}')
    _atomic_csv(final.loc[:, OUTPUT_COLUMNS], output)
    manifest = {
        'reader_schema_version': SCHEMA_VERSION, 'product': PRODUCT,
        'analysis_grade': 'survey_only', 'not_for_science_analysis': True,
        'source': BASE_URL, 'variables': list(args.variables),
        'fields_event_catalog': str(args.events_csv) if args.events_csv else None,
        'fields_event_catalog_sha256': catalog_hash,
        'fields_event_count': event_count, 'merged_event_windows': merged_window_count,
        'before_minutes': args.before_minutes if args.events_csv else None,
        'after_minutes': args.after_minutes if args.events_csv else None,
        'time_selection': 'integration_overlap' if intervals is not None else 'epoch_midpoint',
        'event_windows_csv': str(output.with_suffix('.event_windows.csv')) if event_windows is not None else None,
        'coverage_audit_csv': str(audit_file) if audit_file is not None else None,
        'coverage_audit_sha256': _sha256(audit_file) if audit_file is not None else None,
        'source_days_requested': [d.isoformat() for d in days],
        'source_days_missing': missing, 'source_days_read': source_records,
        'rows_total': int(len(final)),
        'rows_valid': int(final.valid_measurement.astype(str).str.lower().eq('true').sum()),
        'csv_sha256': _sha256(output), 'fingerprint': fingerprint,
        'original_cdfs_deleted': False,
        'warning': 'Summary counts/s only; NOT calibrated differential flux or evidence of particle acceleration.',
    }
    _atomic_text(Path(str(output) + '.manifest.json'), json.dumps(manifest, indent=2) + '\n')
    print(f'[SAVED] {output}: {manifest["rows_total"]} rows, '
          f'{manifest["rows_valid"]} valid; {len(missing)} source days unavailable')
    if missing:
        print('[COVERAGE WARNING] Missing source days are not quiet solar wind:', ', '.join(missing))
    print('[SURVEY ONLY] Do not interpret summary count rates as calibrated particle flux')


def file_mode(args):
    if not args.file.is_file():
        raise FileNotFoundError(args.file)
    frame = read_cdf(args.file, args.variables)
    if args.start:
        frame = frame.loc[pd.to_datetime(frame.timestamp_utc, utc=True) >= _utc(args.start)]
    if args.end:
        frame = frame.loc[pd.to_datetime(frame.timestamp_utc, utc=True) < _utc(args.end)]
    if frame.empty:
        raise RuntimeError('No CDF rate samples in requested UTC interval')
    output = args.output or Path('data/isois/summary_rates_local.csv')
    _atomic_csv(frame, output)
    manifest = {'reader_schema_version': SCHEMA_VERSION, 'product': PRODUCT,
                'analysis_grade': 'survey_only', 'not_for_science_analysis': True,
                'source_files': [args.file.name], 'source_sha256': _sha256(args.file),
                'variables': list(args.variables), 'rows_total': len(frame),
                'rows_valid': int(frame.valid_measurement.sum()),
                'csv_sha256': _sha256(output),
                'warning': 'Survey count rates only, not calibrated flux'}
    _atomic_text(Path(str(output) + '.manifest.json'), json.dumps(manifest, indent=2) + '\n')
    print(f'[SAVED] {output}: {len(frame)} rows')


def date_mode(args):
    if not args.start or not args.end:
        raise ValueError('Supply --start and --end (exclusive), or --events-csv')
    first, last = date.fromisoformat(args.start), date.fromisoformat(args.end)
    if first >= last:
        raise ValueError('--end must be after --start')
    days = [first + timedelta(days=i) for i in range((last - first).days)]
    output = args.output or Path(f'data/isois/summary_rates_{first}_{last}.csv')
    fingerprint = _fingerprint(variables=args.variables, intervals=[],
                               catalog_hash=f'date:{first}:{last}', day_halo=0)
    process_sources(args, days=days, intervals=None, output=output, fingerprint=fingerprint)


def _summary_self_test():
    class FakeCDF:
        def __init__(self, data=None):
            self.data = [5., -1e31, -2., 1e8] if data is None else data
        def varget(self, name):
            return {'A_H_Rate_TS': self.data, 'Epoch': list(range(4)),
                    'Epoch_DELTA': [30e9] * 4}[name]
        def varattsget(self, name):
            return {
                'A_H_Rate_TS': {'DEPEND_0': 'Epoch', 'UNITS': 'counts s!E-1!N',
                                'FILLVAL': -1e31, 'VALIDMIN': 0, 'VALIDMAX': 1e6},
                'Epoch': {'DELTA_PLUS_VAR': 'Epoch_DELTA', 'DELTA_MINUS_VAR': 'Epoch_DELTA'},
                'Epoch_DELTA': {'UNITS': 'ns', 'SI_conversion': '1e-9>seconds'},
            }[name]
    def epoch_converter(_):
        return pd.date_range('2018-11-06', periods=4, freq='min').to_numpy()
    rate = extract_rate(FakeCDF(), 'A_H_Rate_TS',
                        'psp_isois_l2-summary_20181106_v26.cdf', epoch_converter=epoch_converter)
    assert rate.valid_measurement.tolist() == [True, False, False, False]
    assert rate.quality_reason.iloc[1:].tolist() == ['fill_value', 'below_valid_min', 'above_valid_max']
    assert rate.integration_start_utc.iloc[0].startswith('2018-11-05T23:59:30')
    assert rate.integration_end_utc.iloc[0].startswith('2018-11-06T00:00:30')
    try:
        extract_rate(FakeCDF(np.zeros((4, 8))), 'A_H_Rate_TS',
                     'psp_isois_l2-summary_20181106_v26.cdf', epoch_converter=epoch_converter)
    except ValueError as e:
        assert '1D' in str(e)
    else:
        raise AssertionError('Invalid 2D variable silently accepted')
    assert np.isclose(_half_width_seconds([30e9], {'UNITS': 'ns'})[0], 30)
    try:
        _half_width_seconds([30], {'UNITS': 'unknown'})
    except ValueError:
        pass
    else:
        raise AssertionError('Unknown time unit silently accepted')
    events = [('evt-a', _utc('2018-11-06T00:00:00Z'), _utc('2018-11-06T00:00:10Z'), 'test'),
              ('evt-b', _utc('2018-11-06T00:00:20Z'), _utc('2018-11-06T00:00:40Z'), 'test')]
    windows = build_event_windows(events, before_minutes=1, after_minutes=1)
    merged = merged_intervals(windows)
    assert len(merged) == 1
    assert required_days(merged, halo_days=1) == [date(2018, 11, 4), date(2018, 11, 5),
                                                 date(2018, 11, 6), date(2018, 11, 7)]
    # pandas 3.x often parses ISO times as datetime64[us], while Timestamp.value is ns.
    microsecond_times = pd.Series(pd.date_range('2018-11-06', periods=2, freq='min', tz='UTC').as_unit('us'))
    assert _epoch_nanoseconds(microsecond_times).tolist() == [
        _utc('2018-11-06T00:00:00Z').value, _utc('2018-11-06T00:01:00Z').value]
    selected = filter_by_overlap(rate, merged)
    assert len(selected) == 3  # three integrations intersect the merged window
    # Window only touches the integration endpoint: half-open exclusion.
    boundary = [(_utc('2018-11-06T00:00:30Z'), _utc('2018-11-06T00:00:31Z'))]
    assert rate.iloc[0:1].pipe(filter_by_overlap, boundary).empty
    assert source_index('psp_isois_l2-summary_20181106_v4.cdf  '
                        'psp_isois_l2-summary_20181106_v26.cdf')['20181106'].endswith('_v26.cdf')
    # Correct FIELDS 1-second-inclusive -> exclusive conversion and manifest check.
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        catalog = folder / 'mock_strong_events.csv'
        pd.DataFrame([{'event_id': 'e1', 'event_start_utc': '2018-11-06T12:00:00',
                       'event_end_utc': '2018-11-06T12:00:09'}]).to_csv(catalog, index=False)
        manifest = folder / 'mock_strong_retention.json'
        manifest.write_text(json.dumps({'strong_events_retained': 1,
                                        'strong_events_sha256': _sha256(catalog)}))
        e = validate_fields_events(catalog)
        assert e[0][2] == _utc('2018-11-06T12:00:10Z')
        # Invalid checksum is never silently accepted.
        manifest.write_text(json.dumps({'strong_events_retained': 1, 'strong_events_sha256': 'bad'}))
        try:
            validate_fields_events(catalog)
        except ValueError as err:
            assert 'integrity' in str(err)
        else:
            raise AssertionError('Tampered FIELDS catalog was accepted')
    print('PASS: CDF metadata/time, invalid masks, event boundaries and checksum, UTC overlap and halo, version selection')


def _summary_main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--self-test', action='store_true', help='Run synthetic tests; no NASA access')
    parser.add_argument('--inspect', type=Path, help='Print the variable inventory of a local CDF')
    parser.add_argument('--file', type=Path, help='Read one existing summary CDF')
    parser.add_argument('--events-csv', '--fields-events', type=Path,
                        help='FIELDS strong_events or strong_registry CSV (verified against retention JSON)')
    parser.add_argument('--before-minutes', type=float, default=30,
                        help='Context before each magnetic event (default: 30)')
    parser.add_argument('--after-minutes', type=float, default=30,
                        help='Context after each magnetic event (default: 30)')
    parser.add_argument('--day-halo', type=int, choices=(0, 1), default=1,
                        help='Fetch adjacent daily files for integrations straddling UTC midnight')
    parser.add_argument('--allow-unverified-events', action='store_true',
                        help='Explicitly bypass strong-retention checksum (not recommended)')
    parser.add_argument('--plan-only', action='store_true',
                        help='Print event-window plan; no download, processing or output writes')
    parser.add_argument('--start', help='Date-mode UTC start (inclusive), YYYY-MM-DD')
    parser.add_argument('--end', help='Date-mode UTC end (exclusive), YYYY-MM-DD')
    parser.add_argument('--download', action='store_true',
                        help='Allow downloading missing CDFs from the official public archive')
    parser.add_argument('--offline', action='store_true',
                        help='Use cached CDF files only, do not access the archive')
    parser.add_argument('--resume', action='store_true',
                        help='Verify and reuse completed per-day partitions')
    parser.add_argument('--cache', type=Path, default=DEFAULT_CACHE)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--variables', nargs='+', default=list(DEFAULT_VARIABLES))
    args = parser.parse_args()
    if args.self_test:
        _summary_self_test()
        return
    try:
        if args.inspect:
            inspect_file(args.inspect)
        elif args.events_csv:
            if args.file or args.start or args.end:
                parser.error('--events-csv cannot be combined with --file/--start/--end')
            event_mode(args)
        elif args.file:
            if args.plan_only:
                parser.error('--plan-only requires --events-csv')
            file_mode(args)
        else:
            if args.plan_only:
                parser.error('--plan-only requires --events-csv')
            date_mode(args)
    except (ValueError, OSError, RuntimeError, KeyError, ImportError) as exc:
        parser.exit(1, f'[ERROR] {exc}\n')


# === EPI-HI NATIVE-CADENCE CDF READER (previously particle_highres_reader.py) ===
import argparse
import hashlib
import json
import re
from datetime import date, timedelta
from pathlib import Path
from urllib.parse import urljoin

import numpy as np
import pandas as pd


EPIHI_BASE_URL = 'https://spp-isois.sr.unh.edu/data_public/EPIHi/level2/'
EPIHI_PRODUCTS = ('let1-rates10', 'let1-rates60', 'let1-rates300', 'let1-rates3600',
            'let2-rates60', 'het-rates10', 'het-rates60', 'het-rates300',
            'het-rates3600', 'second-rates')
EPIHI_DEFAULT_PRODUCTS = ('let1-rates10', 'let1-rates60', 'het-rates10', 'second-rates')
EPIHI_FILE_RE = re.compile(
    r'^psp_isois-epihi_l2-(' + '|'.join(map(re.escape, EPIHI_PRODUCTS)) +
    r')_(\d{8})_v(\d+)\.cdf$', re.I
)
EPIHI_SCHEMA_VERSION = 2
EPIHI_CHANNEL_GUIDE = {
    'let1-rates10': '10s LET1 H/He energy spectra; He is NOT the original Z6-28 heavy-ion channel',
    'let1-rates60': '60s LET1 species-resolved C, N, O, Ne, Mg, Si, Fe etc.; not identical to the combined Z6-28 summary channel',
    'let1-rates300': '300s LET1 R1/R26 A/B CNO, Ne-Si and Fe-group SECTORED COUNT RATES (not differential flux; not Z6-28 summary)',
    'let1-rates3600': '3600s LET1 A/B proton and heavy-element spectra; coarse background, not resolved switchbacks',
    'let2-rates60': '60s LET2 C-side spectrum; different detector/view, not a LET1 gap replacement',
    'het-rates10': '10s HET H and electron energy spectra plus He',
    'het-rates60': '60s HET rates/energy spectra where published; coarser event measurements',
    'het-rates300': '300s HET A/B CNO, Ne-Si and Fe-group SECTORED COUNT RATES (not differential flux)',
    'het-rates3600': '3600s HET A/B proton/electron spectra; coarse background, not resolved switchbacks',
    'second-rates': '1s LET1/LET2/HET proton/electron count rates and counts; no differential flux or species-resolved energy grid',
}


def _cdf_epoch_ns(raw):
    import cdflib
    t = pd.DatetimeIndex(pd.to_datetime(
        cdflib.cdfepoch.to_datetime(raw), utc=True, errors='coerce'
    )).as_unit('ns')
    if t.isna().any() or t.has_duplicates or not t.is_monotonic_increasing:
        raise ValueError('CDF Epoch has invalid, duplicate, or out-of-order timestamps')
    return t.asi8.copy()


def _bounds(cdf, midpoint):
    attrs = cdf.varattsget('Epoch')
    widths = []
    for key in ('DELTA_MINUS_VAR', 'DELTA_PLUS_VAR'):
        name = _attr(attrs, key)
        if not name:
            raise ValueError(f'Epoch missing {key}; refusing to guess integration times')
        raw = np.asarray(cdf.varget(str(name)))
        if raw.size == 1:
            raw = np.repeat(raw.reshape(-1), len(midpoint))
        if raw.shape != midpoint.shape:
            raise ValueError(f'{name}: invalid integration time shape {raw.shape}')
        widths.append(_half_width_seconds(raw, cdf.varattsget(str(name))))
    start = midpoint - np.rint(widths[0] * 1e9).astype(np.int64)
    end = midpoint + np.rint(widths[1] * 1e9).astype(np.int64)
    if np.any(end <= start):
        raise ValueError('Nonpositive CDF integration interval')
    return start, end


def _overlap(start, end, windows):
    if windows is None:
        return np.ones(len(start), dtype=bool)
    if not windows:
        return np.zeros(len(start), dtype=bool)
    first = np.asarray([a.value for a, _ in windows], dtype=np.int64)
    last = np.asarray([b.value for _, b in windows], dtype=np.int64)
    index = np.searchsorted(last, start, side='right')
    keep = index < len(windows)
    ii = np.flatnonzero(keep)
    keep[ii] = first[index[ii]] < end[ii]
    return keep


def _json(v):
    if isinstance(v, np.ndarray):
        return _json(v.tolist())
    if isinstance(v, (list, tuple)):
        return [_json(x) for x in v]
    if isinstance(v, np.generic):
        return _json(v.item())
    if isinstance(v, bytes):
        return v.decode('utf-8', 'replace')
    return v if isinstance(v, (str, int, float, bool, type(None))) else str(v)


def _safe_array(a):
    out = np.asarray(a)
    if out.dtype.kind in ('S', 'U'):
        return out.astype(str).copy()
    if out.dtype.kind == 'O':
        if not all(isinstance(x, (str, bytes)) for x in out.flat):
            raise ValueError('Unsupported object array in CDF; pickle not permitted')
        return out.astype(str)
    return out.copy()


def _numeric_mask(a, attrs):
    valid = np.isfinite(a)
    fill = _number(_attr(attrs, 'FILLVAL'))
    if fill is not None:
        valid &= ~np.isclose(a, fill, rtol=1e-7, atol=0)
    low, high = _number(_attr(attrs, 'VALIDMIN')), _number(_attr(attrs, 'VALIDMAX'))
    if low is not None:
        valid &= a >= low
    if high is not None:
        valid &= a <= high
    return valid


def decode(cdf, windows=None, epoch_decoder=_cdf_epoch_ns, product=None):
    info = cdf.cdf_info()
    names = list(info.zVariables) + list(info.rVariables)
    if 'Epoch' not in names:
        raise ValueError('Missing CDF Epoch')
    midpoint = epoch_decoder(cdf.varget('Epoch'))
    start, end = _bounds(cdf, midpoint)
    keep = _overlap(start, end, windows)
    arrays = {
        'epoch_utc_ns': midpoint[keep], 'integration_start_utc_ns': start[keep],
        'integration_end_utc_ns': end[keep],
    }
    details, depends, skipped_empty = {}, set(), []
    for name in names:
        if name == 'Epoch':
            continue
        attrs = cdf.varattsget(name)
        if _attr(attrs, 'DEPEND_0') != 'Epoch':
            continue
        a = _safe_array(cdf.varget(name))
        if a.ndim >= 1 and a.shape[0] == 0 and len(midpoint):
            # Some hourly HET CDFs declare empty auxiliary variables,
            # e.g. A_CNO_SECT_Rate with shape (0, 2, 25).
            # Preserve the absence honestly: never pad with zeros.
            skipped_empty.append({'name': name, 'shape': list(a.shape),
                                  'reason': 'CDF variable has zero records'})
            continue
        if a.ndim < 1 or a.shape[0] != len(midpoint):
            raise ValueError(f'{name}: time axis {a.shape} inconsistent with Epoch')
        a = a[keep].copy()
        arrays[name] = a
        # Numeric mask is independent of physical quality flags.
        if np.issubdtype(a.dtype, np.number) and name != 'Quality_Flag' and (
            name.endswith(('_Rate', '_Flux', '_Uncertainty',
                           '_Uncertainty_Plus', '_Uncertainty_Minus'))
            or re.match(r'^(R\d+)?[ABC]_(H|He|Electrons|C|N|O|Fe|Mg|Si)($|_)', name)
            or re.fullmatch(r'(LET[12]_[ABC]|HET_[AB])_(H|Electrons)', name)
        ):
            arrays[name + '__valid_numeric'] = _numeric_mask(a, attrs)
        axes = {key: _attr(attrs, key) for key in
                ('DEPEND_1', 'DEPEND_2', 'DEPEND_3')}
        depends.update(str(value) for value in axes.values() if value)
        details[name] = {
            'shape': list(a.shape), 'units': _json(_attr(attrs, 'UNITS', '')),
            'dependencies': {key: str(value) for key, value in axes.items() if value},
            'fill': _json(_attr(attrs, 'FILLVAL')),
            'valid_min': _json(_attr(attrs, 'VALIDMIN')),
            'valid_max': _json(_attr(attrs, 'VALIDMAX')),
        }
    # Save energy grids, energy widths, sector directions, and quality labels.
    for name in names:
        if name == 'Epoch' or name in arrays:
            continue
        attrs = cdf.varattsget(name)
        if _attr(attrs, 'DEPEND_0'):
            continue
        if name in depends or re.search(r'ENERGY|_LABL|_NUM$|_DELTA', name, re.I):
            a = _safe_array(cdf.varget(name))
            arrays[name] = a
            details[name] = {
                'shape': list(a.shape), 'units': _json(_attr(attrs, 'UNITS', ''))
            }
    flux_names = [n for n in arrays if n.endswith('_Flux')]
    rate_names = [n for n in arrays if n.endswith('_Rate')]
    if product == 'second-rates':
        required = (
            'LET1_A_H_Rate', 'LET1_A_Electrons_Rate',
            'HET_A_H_Rate', 'HET_A_Electrons_Rate',
        )
        missing = [name for name in required if name not in arrays]
        if missing:
            raise ValueError(f'1-second CDF missing mandatory time-aligned rates: {missing}')
        if flux_names:
            raise ValueError('Unexpected differential-flux arrays in second-rates: review schema')
        for name in rate_names:
            units = str(_attr(cdf.varattsget(name), 'UNITS', '')).lower()
            if 'count' not in units or 's!' not in units:
                raise ValueError(f'{name}: expected count rate units, received {units!r}')
        measurement_class = 'counts_and_count_rates_only'
    elif product in ('let1-rates300', 'het-rates300'):
        if flux_names:
            raise ValueError(
                f'{product}: unexpected differential flux variables '
                f'{flux_names}; review product calibration before decoding')
        # Official 300s LET1/HET products provide directional heavy-ion
        # sectored count rates (s^-1), not calibrated differential flux.
        pattern = (
            r'^R(?:1|26)[AB]_(?:CNO|FeGroup|NetoSi)_SECT_Rate$'
            if product == 'let1-rates300' else
            r'^[AB]_(?:CNO|FeGroup|NetoSi)_SECT_Rate$'
        )
        sectored = [name for name in rate_names if re.fullmatch(pattern, name)]
        if not sectored:
            raise ValueError(
                f'{product} has no time-aligned 5-minute heavy-ion sectored '
                f'count-rate arrays (available rate variables: {rate_names}; '
                f'zero-record variables: '
                f'{[x["name"] for x in skipped_empty if x["name"].endswith("_Rate")]})'
            )
        measurement_class = 'sectored_heavy_ion_count_rates_only'
    elif flux_names:
        measurement_class = 'energy_resolved_differential_flux_and_rates'
    else:
        candidates = [(n, str(_attr(cdf.varattsget(n), 'DEPEND_0', '')))
                      for n in names if n.endswith(('_Flux', '_Rate'))]
        raise ValueError(
            f'{product or "science product"} has no decoded flux arrays; '
            f'candidates and DEPEND_0={candidates[:24]}'
        )
    for name, meta in details.items():
        if name not in arrays:
            continue
        for key, target in meta.get('dependencies', {}).items():
            axis = int(key[-1])
            if target not in arrays:
                raise ValueError(f'{name}: missing {key}={target} support axis')
            if arrays[name].ndim <= axis or arrays[target].ndim == 0 \
                    or arrays[name].shape[axis] != len(arrays[target]):
                raise ValueError(f'{name}: {key}={target} dimension mismatch')
    return arrays, {
        'source_time_records': len(midpoint),
        'retained_time_records': int(keep.sum()),
        'median_integration_seconds': float(np.median((end-start)/1e9))
            if len(midpoint) else None,
        'median_cadence_seconds': float(np.median(np.diff(midpoint)/1e9))
            if len(midpoint) > 1 else None,
        'measurement_class': measurement_class,
        'flux_variables': flux_names,
        'rate_variables': rate_names,
        'quality_flags_preserved': 'Quality_Flag' in arrays,
        'quality_flags_interpreted': False,
        'skipped_zero_record_variables': skipped_empty,
        'variables': details,
    }


def _save_npz(path, arrays):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    try:
        with tmp.open('wb') as f:
            np.savez_compressed(f, **arrays)
        tmp.replace(path)
    finally:
        tmp.unlink(missing_ok=True)


def _remote_index(html, product):
    out = {}
    pattern = r'psp_isois-epihi_l2-[\w-]+_\d{8}_v\d+\.cdf'
    for name in set(re.findall(pattern, html, re.I)):
        match = EPIHI_FILE_RE.fullmatch(name)
        if match and match.group(1).lower() == product:
            day, version = match.group(2), int(match.group(3))
            if day not in out or version > out[day][0]:
                out[day] = (version, name)
    return {day: filename for day, (_, filename) in out.items()}


def _download(session, filename, folder, *, attempts=5,
              validator=None, sleeper=None):
    """Retrieve and validate one CDF with retry and HTTP Range partial resume.

    A '.part' file is never treated as a completed CDF. Servers ignoring Range
    requests receive a clean restart; a resumed transfer requires a matching
    Content-Range. Non-retryable HTTP statuses fail promptly.
    """
    import time
    import requests

    if validator is None:
        def validator(source):
            import cdflib
            info = cdflib.CDF(str(source)).cdf_info()
            if not (info.zVariables or info.rVariables):
                raise IOError('Downloaded CDF contains no variables')
    if sleeper is None:
        sleeper = time.sleep

    folder.mkdir(parents=True, exist_ok=True)
    path = folder / filename
    if path.is_file() and path.stat().st_size:
        return path
    part = folder / (filename + '.part')
    url = urljoin(EPIHI_BASE_URL, filename)
    last_error = None

    for attempt in range(1, attempts + 1):
        offset = part.stat().st_size if part.is_file() else 0
        headers = {'Accept-Encoding': 'identity', 'Connection': 'close'}
        if offset:
            headers['Range'] = f'bytes={offset}-'
        try:
            with session.get(url, stream=True, timeout=(20, 90),
                             headers=headers) as response:
                if response.status_code == 416 and offset:
                    # Never promote a .part file on the basis of its CDF header:
                    # partially transferred files can still have valid headers.
                    part.unlink(missing_ok=True)
                    raise IOError('Server rejected partial Range; restarting from zero')

                response.raise_for_status()
                if response.status_code == 206:
                    match = re.fullmatch(
                        r'bytes (\d+)-(\d+)/(\d+|\*)',
                        response.headers.get('Content-Range', ''),
                    )
                    if not offset or not match or int(match.group(1)) != offset:
                        raise IOError('Unsafe HTTP Range response; partial file left intact')
                    total = int(match.group(3)) if match.group(3) != '*' else None
                    mode = 'ab'
                elif response.status_code == 200:
                    # A server may ignore Range: start again, never append a full file.
                    total_raw = response.headers.get('Content-Length', '')
                    total = int(total_raw) if total_raw.isdigit() else None
                    mode = 'wb'
                else:
                    raise IOError(f'Unexpected HTTP status {response.status_code}')

                with part.open(mode) as stream:
                    for chunk in response.iter_content(chunk_size=256 * 1024):
                        if chunk:
                            stream.write(chunk)
            received = part.stat().st_size if part.is_file() else 0
            if received == 0 or (total is not None and received != total):
                raise IOError(
                    f'Incomplete CDF: {received} of {total or "unknown"} bytes'
                )
            validator(part)
            part.replace(path)
            return path
        except requests.HTTPError as exc:
            last_error = exc
            status = exc.response.status_code if exc.response is not None else None
            if status not in (408, 429, 500, 502, 503, 504):
                raise  # e.g. 404: missing product, not a transient outage
        except (requests.ConnectionError, requests.Timeout,
                requests.exceptions.ChunkedEncodingError, IOError, OSError, ValueError) as exc:
            last_error = exc
        if attempt < attempts:
            delay = min(2 ** (attempt - 1), 12)
            print(f'[RETRY DOWNLOAD {attempt}/{attempts}] {filename}: '
                  f'{type(last_error).__name__}; retry in {delay}s', flush=True)
            sleeper(delay)

    raise IOError(
        f'Download failed after {attempts} attempts for {filename}: {last_error}; '
        f'partial file kept at {part}'
    )


def _coverage(windows, starts, ends):
    """Report physical integration overlap, WITHOUT applying instrument quality flags.

    Binary-search relevant integrations instead of scanning all 1-second samples
    for every event; overlapping integrations are counted once in duration.
    """
    starts = np.asarray(starts, dtype=np.int64)
    ends = np.asarray(ends, dtype=np.int64)
    if len(starts) != len(ends):
        raise ValueError('Integration boundary arrays have different lengths')
    order = np.argsort(starts, kind='stable')
    starts, ends = starts[order], ends[order]
    if np.any(ends <= starts):
        raise ValueError('Nonpositive integration interval in coverage audit')
    max_ends = np.maximum.accumulate(ends) if len(ends) else ends
    rows = []
    for w in windows.itertuples(index=False):
        for label, left, right in (
            ('before', w.window_start_utc, w.event_start_utc),
            ('during', w.event_start_utc, w.event_end_exclusive_utc),
            ('after', w.event_end_exclusive_utc, w.window_end_exclusive_utc),
        ):
            a, b = _utc(left).value, _utc(right).value
            if b <= a:
                continue
            lo = int(np.searchsorted(max_ends, a, side='right'))
            hi = int(np.searchsorted(starts, b, side='left'))
            if hi <= lo:
                pieces = []
            else:
                lefts = np.maximum(starts[lo:hi], a)
                rights = np.minimum(ends[lo:hi], b)
                pieces = sorted((int(s), int(e)) for s, e in zip(lefts, rights)
                                if e > s)
            duration, x0, x1 = 0, None, None
            for s, e in pieces:
                if x0 is None:
                    x0, x1 = s, e
                elif s <= x1:
                    x1 = max(x1, e)
                else:
                    duration += x1-x0
                    x0, x1 = s, e
            if x0 is not None:
                duration += x1-x0
            rows.append({
                'event_id': w.event_id, 'segment': label,
                'integration_overlap_seconds': duration / 1e9,
                'coverage_fraction_unfiltered': min(1.0, duration/(b-a)),
                'overlapping_integrations': len(pieces),
            })
    return pd.DataFrame(rows)


def _process(path, product, output_dir, windows, event_df, catalog_hash, resume):
    import cdflib
    match = EPIHI_FILE_RE.fullmatch(path.name)
    if not match or match.group(1).lower() != product:
        raise ValueError(f'Incorrect {product} filename: {path}')
    fingerprint = hashlib.sha256(json.dumps({
        'cdf': _sha256(path), 'catalog': catalog_hash, 'product': product,
        'windows': [[a.isoformat(), b.isoformat()] for a, b in windows]
            if windows is not None else None,
        'schema': EPIHI_SCHEMA_VERSION,
    }, sort_keys=True).encode()).hexdigest()
    # A full-day extraction must not overwrite the subset saved for FIELDS events.
    suffix = f'_events_{fingerprint[:12]}' if event_df is not None else ''
    dest = output_dir / product / (path.stem + suffix + '.npz')
    manifest_file = Path(str(dest) + '.manifest.json')
    if resume and dest.is_file() and manifest_file.is_file():
        manifest = json.loads(manifest_file.read_text())
        coverage = manifest.get('coverage_csv')
        coverage_ok = (not coverage or
            (Path(coverage).is_file() and
             manifest.get('coverage_sha256') == _sha256(coverage)))
        if (manifest.get('fingerprint') == fingerprint
                and manifest.get('output_sha256') == _sha256(dest)
                and coverage_ok):
            return manifest, 'RESUME'
    cdf = cdflib.CDF(str(path))
    arrays, info = decode(cdf, windows, product=product)
    _save_npz(dest, arrays)
    coverage_path = None
    if event_df is not None:
        coverage_path = dest.with_suffix('.coverage.csv')
        _atomic_csv(_coverage(event_df, arrays['integration_start_utc_ns'],
                              arrays['integration_end_utc_ns']), coverage_path)
    manifest = {
        'schema_version': EPIHI_SCHEMA_VERSION, 'product': product,
        'source_cdf': str(path), 'source_sha256': _sha256(path),
        'output_npz': str(dest), 'output_sha256': _sha256(dest),
        'fingerprint': fingerprint, 'events_sha256': catalog_hash,
        'coverage_csv': str(coverage_path) if coverage_path else None,
        'coverage_sha256': _sha256(coverage_path) if coverage_path else None,
        'timestamp_unit': 'UTC Unix nanoseconds',
        'intervals': 'Epoch midpoint +/- DELTA_MINUS/PLUS',
        'channel_mapping': EPIHI_CHANNEL_GUIDE[product],
        'no_rebinning_or_averaging': True,
        'warning': 'Metadata masks do not replace ISʘIS quality-flag interpretation and count statistics.',
        **info,
    }
    _atomic_text(manifest_file, json.dumps(manifest, indent=2) + '\n')
    return manifest, 'READ'


def _epihi_self_test():
    from tempfile import TemporaryDirectory
    from types import SimpleNamespace
    t0 = _utc('2018-11-06T00:00:00Z').value
    ticks = t0 + np.arange(5, dtype=np.int64)*10_000_000_000
    class FakeCDF:
        def __init__(self):
            self.vars = {
                'Epoch': ticks, 'Epoch_DELTA': np.array([5e9]),
                'H_ENERGY': np.array([2., 4.]),
                'Quality_Flag_NUM': np.arange(10),
                'Quality_Flag': np.zeros((5,10), dtype=np.int16),
                'A_H_Flux': np.ones((5,2)),
                'A_H_Rate': np.ones((5,2)),
                'A_H': np.ones((5,2)),
                'A_H_Uncertainty': np.ones((5,2)),
                'HEAVY_ENERGY': np.array([5., 10., 20.]),
                'SECTOR_NUM': np.arange(4),
                'R26A_CNO_SECT_Flux': np.ones((5,3,4)),
            }
        def cdf_info(self):
            return SimpleNamespace(zVariables=list(self.vars), rVariables=[])
        def varget(self, name):
            return self.vars[name]
        def varattsget(self, name):
            if name == 'Epoch':
                return {'DELTA_MINUS_VAR':'Epoch_DELTA','DELTA_PLUS_VAR':'Epoch_DELTA'}
            if name == 'Epoch_DELTA':
                return {'UNITS':'ns'}
            if name.startswith('A_H'):
                return {'DEPEND_0':'Epoch','DEPEND_1':'H_ENERGY','FILLVAL':-1e31}
            if name == 'Quality_Flag':
                return {'DEPEND_0':'Epoch','DEPEND_1':'Quality_Flag_NUM'}
            if name == 'R26A_CNO_SECT_Flux':
                return {'DEPEND_0':'Epoch','DEPEND_1':'HEAVY_ENERGY',
                        'DEPEND_2':'SECTOR_NUM'}
            return {}
    windows = [(_utc('2018-11-06T00:00:05Z'), _utc('2018-11-06T00:00:15Z'))]
    assert _overlap(ticks-5_000_000_000,ticks+5_000_000_000,windows).tolist() == [False,True,False,False,False]
    fake = FakeCDF()
    arrays, meta = decode(fake, windows, epoch_decoder=lambda _: ticks)
    assert arrays['A_H_Flux'].shape == (1,2)
    assert arrays['Quality_Flag'].shape == (1,10)
    assert arrays['R26A_CNO_SECT_Flux'].shape == (1,3,4)
    assert meta['median_integration_seconds'] == 10.0
    assert meta['retained_time_records'] == 1
    fake.vars['A_H_Flux'] = np.ones((5,3))
    try:
        decode(fake, windows, epoch_decoder=lambda _: ticks)
    except ValueError as exc:
        assert 'dimension' in str(exc)
    else:
        raise AssertionError('Missing energy-dimension guard')
    with TemporaryDirectory() as folder:
        path = Path(folder) / 'test.npz'
        _save_npz(path, arrays)
        with np.load(path, allow_pickle=False) as saved:
            assert saved['A_H_Flux'].shape == (1,2)
            assert saved['R26A_CNO_SECT_Flux'].shape == (1,3,4)
    # 1-second coverage must be computed without quadratic Python scanning.
    from types import SimpleNamespace as _SN
    sample_windows = pd.DataFrame([{
        'event_id': 'synthetic', 'window_start_utc': '2018-11-06T00:00:00Z',
        'event_start_utc': '2018-11-06T00:00:01Z',
        'event_end_exclusive_utc': '2018-11-06T00:00:03Z',
        'window_end_exclusive_utc': '2018-11-06T00:00:04Z'
    }])
    overlap_ns = np.array([t0, t0+1_000_000_000,
                           t0+2_000_000_000, t0+3_000_000_000])
    audit = _coverage(sample_windows, overlap_ns, overlap_ns+1_000_000_000)
    assert audit['segment'].tolist() == ['before', 'during', 'after']
    assert audit['integration_overlap_seconds'].tolist() == [1.0, 2.0, 1.0]
    assert (audit['coverage_fraction_unfiltered'] == 1.0).all()

    class OneSecondCDF(FakeCDF):
        def __init__(self):
            super().__init__()
            t = t0 + np.arange(5, dtype=np.int64)*1_000_000_000
            self.vars = {
                'Epoch': t, 'Epoch_DELTA': np.array([500_000_000], dtype=np.int64),
                'Quality_Flag_NUM': np.arange(10),
                'Quality_Flag': np.zeros((5, 10), dtype=np.int16),
            }
            for instrument, side in (('LET1','A'),('LET1','B'),('LET2','C'),
                                     ('HET','A'),('HET','B')):
                for species in ('H','Electrons'):
                    name = f'{instrument}_{side}_{species}'
                    self.vars[name] = np.ones(5, dtype=np.float64)
                    self.vars[name + '_Rate'] = np.ones(5, dtype=np.float64)
            self.vars['LET1_A_H_Rate'][1] = -1e31

        def varattsget(self, name):
            if name.startswith(('LET1_', 'LET2_', 'HET_')):
                return {
                    'DEPEND_0': 'Epoch',
                    'UNITS': 'counts s!E-1!N' if name.endswith('_Rate') else 'counts',
                    'FILLVAL': -1e31, 'VALIDMIN': 0,
                }
            return super().varattsget(name)

    onesec = OneSecondCDF()
    quicklook, quick_meta = decode(
        onesec, epoch_decoder=lambda _: onesec.vars['Epoch'], product='second-rates'
    )
    assert quick_meta['measurement_class'] == 'counts_and_count_rates_only'
    assert quick_meta['median_integration_seconds'] == 1.0
    assert len(quick_meta['rate_variables']) == 10
    assert quicklook['LET1_A_H_Rate'].shape == (5,)
    assert quicklook['LET1_A_H_Rate__valid_numeric'].tolist() == [
        True, False, True, True, True]
    assert quicklook['LET1_A_H__valid_numeric'].all()
    assert quicklook['Quality_Flag'].shape == (5, 10)
    with TemporaryDirectory() as folder:
        target = Path(folder) / 'second-rates.npz'
        _save_npz(target, quicklook)
        with np.load(target, allow_pickle=False) as saved:
            assert saved['HET_A_Electrons_Rate'].shape == (5,)
    try:
        decode(onesec, epoch_decoder=lambda _: onesec.vars['Epoch'],
               product='het-rates10')
    except ValueError as exc:
        assert 'no decoded flux' in str(exc)
    else:
        raise AssertionError('Rates-only quicklook silently passed as differential flux')
    # Five-minute LET1/HET CDFs have directional heavy-ion count rates.
    # Verify axes, fill masks, cadence and flux/rate discrimination.
    class FiveMinuteCDF(FakeCDF):
        def __init__(self, kind):
            super().__init__()
            self.vars = {
                'Epoch': t0 + np.arange(5, dtype=np.int64)*300_000_000_000,
                'Epoch_DELTA': np.array([150e9]),
                'GROUP_ENERGY': np.array([8., 16.]),
                'GROUP_SECTORS': np.arange(25),
                'Quality_Flag_NUM': np.arange(10),
                'Quality_Flag': np.zeros((5, 10), dtype=np.int16),
            }
            self.rate_name = ('R26A_CNO_SECT_Rate'
                              if kind == 'let1-rates300'
                              else 'A_CNO_SECT_Rate')
            self.vars[self.rate_name] = np.ones((5, 2, 25))
            self.vars[self.rate_name][2, 1, 3] = -1e31
        def varattsget(self, name):
            if name == self.rate_name:
                return {
                    'DEPEND_0': 'Epoch', 'DEPEND_1': 'GROUP_ENERGY',
                    'DEPEND_2': 'GROUP_SECTORS', 'UNITS': 's!E-1!N',
                    'FILLVAL': -1e31,
                }
            return super().varattsget(name)

    for five_product in ('let1-rates300', 'het-rates300'):
        fake_five = FiveMinuteCDF(five_product)
        observed, five_meta = decode(
            fake_five, epoch_decoder=lambda _: fake_five.vars['Epoch'],
            product=five_product)
        assert five_meta['measurement_class'] == 'sectored_heavy_ion_count_rates_only'
        assert five_meta['flux_variables'] == []
        assert five_meta['median_integration_seconds'] == 300.0
        assert fake_five.rate_name in five_meta['rate_variables']
        assert observed[fake_five.rate_name].shape == (5, 2, 25)
        assert not observed[fake_five.rate_name+'__valid_numeric'][2, 1, 3]
        assert observed[fake_five.rate_name+'__valid_numeric'][0, 0, 0]
        assert observed['GROUP_ENERGY'].shape == (2,)
        assert observed['GROUP_SECTORS'].shape == (25,)
        try:
            decode(fake_five, epoch_decoder=lambda _: fake_five.vars['Epoch'],
                   product='het-rates10')
        except ValueError as exc:
            assert 'no decoded flux' in str(exc)
        else:
            raise AssertionError('Five-minute rate mislabeled as 10s flux')

    # Official hourly HET files can declare zero-record auxiliaries.
    # They must be recorded as unavailable, not replaced with fake zeros.
    class HourlyEmptyCDF(FakeCDF):
        def __init__(self):
            super().__init__()
            self.vars['A_CNO_SECT_Rate'] = np.empty((0, 2, 25))
        def varattsget(self, name):
            if name == 'A_CNO_SECT_Rate':
                return {'DEPEND_0':'Epoch','UNITS':'counts s!E-1!N'}
            return super().varattsget(name)
    hourly = HourlyEmptyCDF()
    result, meta = decode(hourly, epoch_decoder=lambda _: ticks,
                          product='het-rates3600')
    assert 'A_CNO_SECT_Rate' not in result
    assert any(x['name'] == 'A_CNO_SECT_Rate'
               for x in meta['skipped_zero_record_variables'])
    assert 'A_H_Flux' in result

    # Read-only HET CDF inspection must distinguish declared, zero-record
    # sectored rate arrays from actual observational coverage.
    with TemporaryDirectory() as folder:
        from contextlib import redirect_stdout
        from unittest.mock import patch
        import io
        target = Path(folder) / 'het-rates300-test.cdf'
        target.write_bytes(b'fake CDF for injected reader')
        output = io.StringIO()
        with patch('cdflib.CDF', return_value=hourly):
            with redirect_stdout(output):
                assert _epihi_main(['--inspect-cdf', str(target)]) == 0
        assert 'A_CNO_SECT_Rate' in output.getvalue()
        assert 'ZERO RECORDS' in output.getvalue()
        assert 'Epoch records: 5' in output.getvalue()
        assert list(Path(folder).iterdir()) == [target]

    # Synthetic transient disconnect and safe HTTP Range resume, no live network.
    import requests
    class FakeResponse:
        def __init__(self, status, payload, headers):
            self.status_code = status
            self.payload = payload
            self.headers = headers
        def __enter__(self):
            return self
        def __exit__(self, exc_type, exc, tb):
            return False
        def raise_for_status(self):
            if self.status_code >= 400:
                raise requests.HTTPError(f'HTTP {self.status_code}')
        def iter_content(self, chunk_size):
            yield self.payload

    class FakeSession:
        def __init__(self, sequence):
            self.sequence = iter(sequence)
            self.calls = []
        def get(self, url, **kwargs):
            self.calls.append(kwargs)
            item = next(self.sequence)
            if isinstance(item, Exception):
                raise item
            return item

    with TemporaryDirectory() as folder:
        folder = Path(folder)
        sample = 'psp_isois-epihi_l2-het-rates10_20181106_v26.cdf'
        part = folder / (sample + '.part')
        part.write_bytes(b'ab')
        session = FakeSession([
            requests.ConnectionError('simulated network disconnect'),
            FakeResponse(206, b'cde', {'Content-Range': 'bytes 2-4/5'}),
        ])
        def verify_sample(path):
            assert path.read_bytes() == b'abcde'
        result = _download(session, sample, folder, attempts=2,
                           validator=verify_sample, sleeper=lambda _: None)
        assert result.read_bytes() == b'abcde'
        assert session.calls[0]['headers']['Range'] == 'bytes=2-'
        assert not part.exists()
        result.unlink()
        part.write_bytes(b'old')
        session = FakeSession([
            FakeResponse(200, b'fresh', {'Content-Length': '5'})
        ])
        result = _download(session, sample, folder, attempts=1,
                           validator=lambda p: None, sleeper=lambda _: None)
        assert result.read_bytes() == b'fresh'  # ignores Range: clean restart

    assert _remote_index(
        'psp_isois-epihi_l2-let1-rates3600_20181206_v26.cdf',
        'let1-rates3600')['20181206'].endswith('_v26.cdf')
    assert _remote_index(
        'psp_isois-epihi_l2-het-rates3600_20181206_v26.cdf',
        'het-rates3600')['20181206'].endswith('_v26.cdf')
    assert EPIHI_FILE_RE.fullmatch(
        'psp_isois-epihi_l2-let1-rates3600_20181206_v26.cdf')
    assert _remote_index('psp_isois-epihi_l2-het-rates10_20181106_v2.cdf '
                         'psp_isois-epihi_l2-het-rates10_20181106_v26.cdf',
                         'het-rates10')['20181106'].endswith('_v26.cdf')
    print('PASS: EPI-Hi timing, energy axes, 1s counts, flags, NPZ, network retry and safe Range resume')


def inspect_epihi_cdf(path):
    """Read-only inventory for a local EPI-Hi CDF; never invent data rows."""
    import cdflib
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f'CDF not found: {path}')
    cdf = cdflib.CDF(str(path))
    info = cdf.cdf_info()
    names = list(info.zVariables) + list(info.rVariables)
    print(f'[CDF INSPECT] {path}')
    if 'Epoch' not in names:
        raise ValueError('CDF lacks Epoch time coordinate')
    epoch = np.asarray(cdf.varget('Epoch'))
    print(f'  Epoch records: {len(epoch)}')
    found = 0
    for name in names:
        if not name.endswith(('_Rate', '_Flux')):
            continue
        attrs = cdf.varattsget(name)
        dependency = _attr(attrs, 'DEPEND_0')
        value = np.asarray(cdf.varget(name))
        records = int(value.shape[0]) if value.ndim else 0
        state = ('ZERO RECORDS' if records == 0 else
                 'TIMING MISMATCH' if dependency == 'Epoch' and records != len(epoch)
                 else 'time-aligned' if dependency == 'Epoch' else 'non-Epoch')
        print(f'  {name}: shape={tuple(value.shape)}; '
              f'DEPEND_0={dependency}; {state}')
        found += 1
    print(f'[CDF INSPECT] {found} rate/flux variables; '
          'inspection only, no data created')
    return 0


def _epihi_main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--self-test', action='store_true')
    p.add_argument('--file', type=Path, help='Decode one existing EPI-Hi CDF')
    p.add_argument('--inspect-cdf', type=Path,
                   help='Read-only inventory of Epoch and real rate/flux array shapes')
    p.add_argument('--events-csv', type=Path)
    p.add_argument('--products', nargs='+', choices=EPIHI_PRODUCTS, default=list(EPIHI_DEFAULT_PRODUCTS),
                   help='Select EPI-Hi products (1s rates-only quicklook; 10s/60s energy-resolved science)')
    p.add_argument('--start', help='UTC YYYY-MM-DD inclusive')
    p.add_argument('--end', help='UTC YYYY-MM-DD exclusive')
    p.add_argument('--before-seconds', type=float, default=10,
                   help='Immediate pre-event observation window (default: 10 s)')
    p.add_argument('--after-seconds', type=float, default=10,
                   help='Immediate post-event observation window (default: 10 s)')
    p.add_argument('--before-minutes', type=float, default=None,
                   help='Legacy option overriding --before-seconds')
    p.add_argument('--after-minutes', type=float, default=None,
                   help='Legacy option overriding --after-seconds')
    p.add_argument('--day-halo', type=int, choices=(0,1), default=1)
    p.add_argument('--cache', type=Path, default=Path('data/isois/epihi_cdf'))
    p.add_argument('--output-dir', type=Path, default=Path('data/isois/highres'))
    p.add_argument('--download', action='store_true')
    p.add_argument('--offline', action='store_true')
    p.add_argument('--resume', action='store_true')
    p.add_argument('--plan-only', action='store_true')
    p.add_argument('--download-attempts', type=int, default=4,
                   help='Retries per source CDF for transient network errors (default: 4)')
    args = p.parse_args(argv)
    if args.inspect_cdf:
        return inspect_epihi_cdf(args.inspect_cdf)
    if args.self_test:
        _epihi_self_test()
        return
    if args.offline and args.download:
        p.error('--offline and --download are mutually exclusive')
    if not 1 <= args.download_attempts <= 12:
        p.error('--download-attempts must be between 1 and 12')
    before_seconds = (args.before_minutes*60 if args.before_minutes is not None
                      else args.before_seconds)
    after_seconds = (args.after_minutes*60 if args.after_minutes is not None
                     else args.after_seconds)
    if not np.isfinite(before_seconds) or not np.isfinite(after_seconds) \
            or min(before_seconds, after_seconds) < 0:
        p.error('Invalid event padding')
    event_df, windows, catalog_hash = None, None, None
    if args.events_csv:
        if args.start or args.end:
            p.error('--events-csv cannot be combined with --start/--end')
        events = validate_fields_events(args.events_csv)
        event_df = build_event_windows(
            events, before_minutes=before_seconds/60, after_minutes=after_seconds/60
        )
        print(f'[EVENT WINDOWS] {len(events)} candidates, '
              f'{before_seconds:g}s before + {after_seconds:g}s after each event',
              flush=True)
        windows = merged_intervals(event_df)
        catalog_hash = _sha256(args.events_csv)
    results, missing, failures = [], [], []

    def checkpoint(*, complete=False):
        """Create an explicit partial manifest after every completed or failed day."""
        _atomic_text(args.output_dir / 'index.json', json.dumps({
            'schema_version': EPIHI_SCHEMA_VERSION,
            'status': 'complete' if complete else 'partial',
            'fields_catalog': str(args.events_csv) if args.events_csv else None,
            'fields_catalog_sha256': catalog_hash,
            'event_padding_seconds': {'before': before_seconds, 'after': after_seconds},
            'products': ([EPIHI_FILE_RE.fullmatch(args.file.name).group(1)]
                         if args.file else args.products),
            'files': results,
            'missing_source_product_days': missing,
            'failed_source_product_days': failures,
            'warning': ('Time coverage does not imply good particle statistics, '
                        'valid quality flags, or association with switchbacks. '
                        'One-second counts/rates are not energy-resolved flux.'),
        }, indent=2) + '\n')

    def process(product, path):
        try:
            manifest, action = _process(
                path, product, args.output_dir, windows,
                event_df, catalog_hash, args.resume
            )
        except (ValueError, OSError, KeyError, RuntimeError) as exc:
            record = {'product': product, 'source_file': str(path),
                      'error': f'{type(exc).__name__}: {exc}'}
            failures.append(record)
            print(f'[PROCESS FAILED] {product} {path.name}: '
                  f'{record["error"]}', flush=True)
            checkpoint()
            return
        print(f'[{action} {len(results)+1}] {product}: '
              f'{manifest["retained_time_records"]} integrations; '
              f'{manifest["median_integration_seconds"]}s duration',
              flush=True)
        results.append({
            'product': product, 'source_file': str(path),
            'npz': manifest['output_npz'],
            'time_records': manifest['retained_time_records'],
        })
        checkpoint()

    if args.file:
        if args.plan_only:
            p.error('--plan-only cannot be used with --file')
        if not args.file.is_file():
            p.error(f'Local file does not exist: {args.file}')
        match = EPIHI_FILE_RE.fullmatch(args.file.name)
        if not match:
            p.error('Unsupported high-resolution EPI-Hi filename')
        process(match.group(1).lower(), args.file)
    else:
        if windows is not None:
            days = required_days(windows, halo_days=args.day_halo)
        else:
            if not args.start or not args.end:
                p.error('Provide --events-csv, --file, or --start/--end')
            first, last = date.fromisoformat(args.start), date.fromisoformat(args.end)
            if last <= first:
                p.error('--end must be after --start')
            days = [first + timedelta(days=i) for i in range((last-first).days)]
        print(f'[HIGHRES PLAN] {len(windows) if windows is not None else 0} merged windows; '
              f'{len(days)} UTC source days; products: {", ".join(args.products)}',
              flush=True)
        if args.plan_only:
            print('[PLAN ONLY] No remote access or output writes')
            return

        import requests
        import time
        with requests.Session() as session:
            remote = {}
            archive_error = None
            if not args.offline:
                for attempt in range(1, args.download_attempts+1):
                    try:
                        response = session.get(
                            EPIHI_BASE_URL, timeout=(20, 60),
                            headers={'Connection': 'close'}
                        )
                        response.raise_for_status()
                        remote = {name: _remote_index(response.text, name)
                                  for name in args.products}
                        break
                    except requests.RequestException as exc:
                        archive_error = str(exc)
                        if attempt < args.download_attempts:
                            delay = min(2 ** (attempt-1), 12)
                            print(f'[RETRY ARCHIVE {attempt}/{args.download_attempts}] '
                                  f'{type(exc).__name__}; retry in {delay}s',
                                  flush=True)
                            time.sleep(delay)
                if not remote:
                    print('[ARCHIVE INDEX UNAVAILABLE] '
                          f'{archive_error}; attempting cached files only',
                          flush=True)

            # Process each daily file as soon as it is acquired, so a later
            # interrupted download never discards completed partitions.
            for product in args.products:
                folder = args.cache / product
                local = {}
                for cached in folder.glob('*.cdf'):
                    m = EPIHI_FILE_RE.fullmatch(cached.name)
                    if (m and m.group(1).lower() == product
                            and cached.stat().st_size):
                        key, version = m.group(2), int(m.group(3))
                        if key not in local or version > local[key][0]:
                            local[key] = (version, cached)

                for day in days:
                    key = day.strftime('%Y%m%d')
                    name = remote.get(product, {}).get(key)
                    remote_path = folder / name if name else None
                    cached_path = local.get(key, (None, None))[1]
                    # If index is available, its newest version wins. Offline
                    # or missing-index days may use the verified local cache.
                    path = (remote_path if remote_path is not None
                            else cached_path)
                    if path is None:
                        missing.append(f'{product} {day}: no source CDF listed or cached')
                        checkpoint()
                        continue
                    if not path.exists():
                        if not args.download:
                            missing.append(f'{product} {day}: not cached; use --download')
                            checkpoint()
                            continue
                        try:
                            path = _download(
                                session, path.name, folder,
                                attempts=args.download_attempts,
                            )
                        except (requests.RequestException, OSError) as exc:
                            failures.append({
                                'product': product, 'source_file': path.name,
                                'error': f'{type(exc).__name__}: {exc}',
                            })
                            print(f'[DOWNLOAD FAILED] {product} {day}: {exc}',
                                  flush=True)
                            checkpoint()
                            continue
                    process(product, path)

    if event_df is not None:
        _atomic_csv(event_df, args.output_dir / 'fields_event_windows.csv')
    complete = bool(results) and not failures and not missing
    checkpoint(complete=complete)
    print(f'[SAVED] {len(results)} native-cadence partitions; '
          f'{len(missing)} missing; {len(failures)} failed; '
          f'index: {args.output_dir / "index.json"}', flush=True)
    if not complete:
        print('[INCOMPLETE] Event windows with missing source products must not '
              'be interpreted as quiet solar wind. Re-run with --resume '
              'to retry remaining product-days.', flush=True)
        raise SystemExit(2)


# ---------------------------------------------------------------------------
# Unified public API (data ingestion only, no particle science calculations)
# ---------------------------------------------------------------------------
def list_products():
    print("NASA IS⊙IS L2 SUMMARY: separate hourly survey rates")
    print("NASA EPI-Hi native-cadence products:")
    for product in EPIHI_PRODUCTS:
        default = " [default]" if product in EPIHI_DEFAULT_PRODUCTS else ""
        print(f"  {product}{default}: {EPIHI_CHANNEL_GUIDE[product]}")
    return 0


def self_test():
    _summary_self_test()
    _epihi_self_test()
    # Public dispatch must never implicitly set --download.
    assert "--download" not in ["--self-test"]
    assert "het-rates3600" in EPIHI_PRODUCTS
    print("PASS: both IS⊙IS reader decoders and unified ingestion interface")


def help_text():
    print("""Unified IS⊙IS reader — real CDF decoding in this single source file

  python project.py isois-reader epihi [EPI-Hi options]
  python project.py isois-reader summary [summary survey options]
  python project.py isois-reader list-products
  python project.py isois-reader --self-test

Use 'epihi --help' or 'summary --help' for full compatible options.
Use 'epihi --inspect-cdf /path/to/local.cdf' to inspect variable record
counts and shapes read-only (useful for empty HET 5-minute products).
All cached source CDFs, native NPZs and summary CSV formats stay unchanged.
No download happens unless --download is explicitly requested.
""")


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("--help", "-h", "help"):
        help_text()
        return
    if argv == ["--self-test"]:
        self_test()
        return
    if argv[0] in ("list-products", "--list-products"):
        if len(argv) != 1:
            raise SystemExit("list-products takes no arguments")
        list_products()
        return
    if argv[0] in ("--highres", "highres"):
        argv = ["epihi", *argv[1:]]
    mode, *args = argv
    if mode == "epihi":
        return _epihi_main(args)
    if mode == "summary":
        # Summary's legacy parser reads sys.argv, unlike EPI-Hi's argv parser.
        saved = sys.argv
        try:
            sys.argv = [saved[0], *args]
            return _summary_main()
        finally:
            sys.argv = saved
    raise SystemExit(f"Unknown reader mode {mode!r}; use epihi or summary")


if __name__ == "__main__":
    main()
