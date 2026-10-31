#!/usr/bin/env python3
"""Parker Solar Probe ISʘIS L2 *summary-rate survey* reader.

Default: read exact detector integration intervals from a daily summary CDF.
Event mode: use existing FIELDS strong-event UTC times to retrieve only needed
ISʘIS days, retain every particle integration overlapping a padded event, and
write an audit and reusable event-window catalog. No FIELDS reanalysis occurs.

IMPORTANT: L2 summary count rates are survey data, NOT calibrated differential
flux, energy spectra, pitch angles or a science-grade correlation result.
No raw CDFs or previously generated FIELDS files are ever deleted.
Run from the repository root so relative data paths resolve correctly.
"""
from __future__ import annotations

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
        raise RuntimeError('Missing cdflib: python -m pip install -r code/M_FIELD/external_libraries_requirements.txt') from exc
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
    data = pd.to_numeric(pd.Series(values), errors='coerce').to_numpy(dtype=float)
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
    left = pd.to_datetime(frame['integration_start_utc'], utc=True).astype('int64').to_numpy()
    right = pd.to_datetime(frame['integration_end_utc'], utc=True).astype('int64').to_numpy()
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
        group['_start_ns'] = pd.to_datetime(group.integration_start_utc, utc=True).astype('int64')
        group['_end_ns'] = pd.to_datetime(group.integration_end_utc, utc=True).astype('int64')
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


def self_test():
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


def main():
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
        self_test()
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


if __name__ == '__main__':
    main()
