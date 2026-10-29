#!/usr/bin/env python3
"""PSP ISʘIS *summary-product* survey reader (NOT science-grade flux).

Reads actual CDF metadata: each rate's DEPEND_0 epoch, DELTA_PLUS/MINUS_VAR
integration half widths, FILLVAL, VALIDMIN/MAX and UNITS. Missing data remain
missing. No assumption that a timestamp represents an instantaneous measurement.

The ISʘIS team's summary-product warning is important: use summary rates for
survey/triage, then use the appropriate EPI-Lo/EPI-Hi science products for
physical flux, spectra, pitch angle or publication-grade inference.
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
    'A_H_Rate_TS',
    'A_Heavy_Rate_TS',
    'HET_A_H_Rate_TS',
    'HET_A_Electrons_Rate_TS',
)
DEFAULT_CACHE = Path('data/isois/summary_cdf')
OUTPUT_COLUMNS = [
    'timestamp_utc', 'integration_start_utc', 'integration_end_utc',
    'variable', 'value', 'units', 'valid_measurement', 'quality_reason',
    'source_file', 'source_version', 'product', 'analysis_grade',
]


def _cdf_module():
    try:
        import cdflib
    except ImportError as exc:
        raise RuntimeError('Missing cdflib; install via: python -m pip install cdflib') from exc
    return cdflib


def _attr(attrs, name, default=None):
    value = attrs.get(name, default)
    if isinstance(value, np.ndarray) and value.size == 1:
        return value.reshape(-1)[0].item()
    if isinstance(value, (tuple, list)) and len(value) == 1:
        return value[0]
    return value


def _number(value):
    if value is None:
        return None
    try:
        z = float(np.asarray(value).reshape(-1)[0])
    except (ValueError, TypeError, IndexError):
        return None
    return z if np.isfinite(z) else None


def _iso(ts):
    return ts.isoformat().replace('+00:00', 'Z')


def _utc_times(raw, converter):
    converted = converter(raw)
    times = pd.to_datetime(converted, utc=True, errors='coerce')
    return pd.DatetimeIndex(times)


def _half_width_seconds(delta, attrs):
    """Return the actual half-integration duration, not guessed cadence."""
    units = str(_attr(attrs, 'UNITS', '')).lower().strip()
    si = str(_attr(attrs, 'SI_conversion', '')).lower().strip()
    if 'nanosecond' in units or units in ('ns', 'nsec') or '1e-9>seconds' in si:
        scale = 1e-9
    elif 'millisecond' in units or units in ('ms', 'msec') or '1e-3>seconds' in si:
        scale = 1e-3
    elif units in ('s', 'sec', 'second', 'seconds') or '1>seconds' in si:
        scale = 1.0
    else:
        raise ValueError(f'Unknown epoch delta units={units!r}, SI_conversion={si!r}')
    result = np.asarray(delta, dtype=float) * scale
    if result.ndim != 1 or np.any(~np.isfinite(result)) or np.any(result <= 0):
        raise ValueError('Invalid positive epoch half-width values')
    return result


def extract_rate(cdf, variable, source_file, *, epoch_converter):
    """Extract a 1D rate variable with truthful per-sample integration times.

    Accepts an injected CDF-like object for deterministic synthetic tests.
    Rejects unlabelled/2D variables instead of silently collapsing angles or
    energy bins. Use product-specific loaders for spectra/pitch-angle data.
    """
    attrs = cdf.varattsget(variable)
    epoch_name = _attr(attrs, 'DEPEND_0')
    units = str(_attr(attrs, 'UNITS', '')).strip()
    if not epoch_name or not units:
        raise ValueError(f'{variable}: missing DEPEND_0 or UNITS')
    if 'rate' not in variable.lower() or ('count' not in units.lower() and 's!' not in units.lower()):
        raise ValueError(f'{variable}: only explicitly labelled count rates allowed; units={units!r}')
    values = np.asarray(cdf.varget(variable))
    if values.ndim != 1:
        raise ValueError(f'{variable}: expected 1D rates; got shape {values.shape}. Do not flatten energy/look bins.')
    epochs = _utc_times(cdf.varget(str(epoch_name)), epoch_converter)
    if len(values) != len(epochs):
        raise ValueError(f'{variable}: epoch/data length mismatch')
    epoch_attrs = cdf.varattsget(str(epoch_name))
    plus_name = _attr(epoch_attrs, 'DELTA_PLUS_VAR')
    minus_name = _attr(epoch_attrs, 'DELTA_MINUS_VAR')
    if not plus_name or not minus_name:
        raise ValueError(f'{variable}: missing integration delta metadata on {epoch_name}')

    def width(name):
        raw = np.asarray(cdf.varget(str(name)))
        if raw.ndim == 0:
            raw = np.repeat(raw, len(values))
        if raw.shape != (len(values),):
            raise ValueError(f'{name}: expected one integration half-width per sample')
        return _half_width_seconds(raw, cdf.varattsget(str(name)))

    plus, minus = width(plus_name), width(minus_name)
    start = epochs - pd.to_timedelta(minus, unit='s')
    end = epochs + pd.to_timedelta(plus, unit='s')
    if epochs.isna().any() or start.isna().any() or end.isna().any():
        raise ValueError(f'{variable}: unparseable epoch time')
    if not epochs.is_monotonic_increasing or epochs.has_duplicates:
        raise ValueError(f'{variable}: epochs not strictly increasing')
    vals = pd.to_numeric(pd.Series(values), errors='coerce').to_numpy(dtype=float)
    invalid = ~np.isfinite(vals)
    reasons = np.full(len(vals), '', dtype=object)
    reasons[invalid] = 'nonfinite'
    fill = _number(_attr(attrs, 'FILLVAL'))
    if fill is not None:
        marked = np.isclose(vals, fill, rtol=1e-7, atol=0, equal_nan=False)
        invalid |= marked
        reasons[marked] = 'fill_value'
    lower = _number(_attr(attrs, 'VALIDMIN'))
    upper = _number(_attr(attrs, 'VALIDMAX'))
    if lower is not None:
        marked = np.isfinite(vals) & (vals < lower) & ~invalid
        invalid |= marked
        reasons[marked] = 'below_valid_min'
    if upper is not None:
        marked = np.isfinite(vals) & (vals > upper) & ~invalid
        invalid |= marked
        reasons[marked] = 'above_valid_max'
    marked = np.isfinite(vals) & (vals < 0) & ~invalid
    invalid |= marked
    reasons[marked] = 'negative_rate'
    vals[invalid] = np.nan
    match = FILE_PATTERN.fullmatch(Path(source_file).name)
    if not match:
        raise ValueError(f'Unexpected ISʘIS summary filename: {Path(source_file).name}')
    return pd.DataFrame({
        'timestamp_utc': [_iso(x) for x in epochs],
        'integration_start_utc': [_iso(x) for x in start],
        'integration_end_utc': [_iso(x) for x in end],
        'variable': variable, 'value': vals, 'units': units,
        'valid_measurement': ~invalid, 'quality_reason': reasons,
        'source_file': Path(source_file).name, 'source_version': int(match.group(2)),
        'product': PRODUCT, 'analysis_grade': 'survey_only',
    }, columns=OUTPUT_COLUMNS)


def inspect_file(path):
    """Print CDF variable inventory, dimensions and essential ISTP attributes."""
    cdf = _cdf_module().CDF(str(path))
    info = cdf.cdf_info()
    variables = list(info.zVariables) + list(info.rVariables)
    for name in variables:
        attrs = cdf.varattsget(name)
        description = str(_attr(attrs, 'CATDESC', ''))
        print(f'{name:45s} depend={str(_attr(attrs, "DEPEND_0", "")):18s} '
              f'units={str(_attr(attrs, "UNITS", "")):22s} {description[:75]}')
    print(f'[INVENTORY] {len(variables)} CDF variables in {path}')


def read_cdf(path, variables):
    cdflib = _cdf_module()
    cdf = cdflib.CDF(str(path))
    available_info = cdf.cdf_info()
    available = set(available_info.zVariables) | set(available_info.rVariables)
    missing = set(variables) - available
    if missing:
        raise ValueError(f'{path.name}: requested variables absent: {sorted(missing)}')
    return pd.concat([
        extract_rate(cdf, name, path.name, epoch_converter=cdflib.cdfepoch.to_datetime)
        for name in variables
    ], ignore_index=True)


def _source_index(session):
    response = session.get(BASE_URL, timeout=(15, 90))
    response.raise_for_status()
    found = {}
    for day, version in re.findall(r'psp_isois_l2-summary_(\d{8})_v(\d+)\.cdf', response.text):
        filename = f'psp_isois_l2-summary_{day}_v{version}.cdf'
        old = found.get(day)
        if old is None or int(version) > old[0]:
            found[day] = (int(version), filename)
    if not found:
        raise RuntimeError('No expected summary CDF filenames found in ISʘIS source listing')
    return {day: filename for day, (_, filename) in found.items()}


def _download(session, filename, cache, *, attempts=4):
    cache.mkdir(parents=True, exist_ok=True)
    dest = cache / filename
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    part = dest.with_name(dest.name + '.part')
    for trial in range(1, attempts + 1):
        part.unlink(missing_ok=True)
        try:
            with session.get(urljoin(BASE_URL, filename), timeout=(20, 120), stream=True) as response:
                response.raise_for_status()
                expected = int(response.headers.get('Content-Length', 0))
                total = 0
                with part.open('wb') as target:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            target.write(chunk)
                            total += len(chunk)
                if total == 0 or (expected > 0 and total != expected):
                    raise IOError(f'Incomplete transfer: got {total} bytes, expected {expected}')
            part.replace(dest)
            return dest
        except Exception as exc:
            print(f'[RETRY {trial}/{attempts}] {filename}: {exc}', file=sys.stderr)
            if trial == attempts:
                raise
            time.sleep(min(trial * 2, 8))
        finally:
            part.unlink(missing_ok=True)
    raise AssertionError('unreachable')


def _digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def _parse_day(value):
    return date.fromisoformat(value)


def run(args):
    if args.inspect:
        inspect_file(args.inspect)
        return
    if args.file:
        paths = [args.file]
        interval_start = _parse_day(args.start) if args.start else None
        interval_end = _parse_day(args.end) if args.end else None
        missing_days = []
    else:
        if not args.start or not args.end or not args.download:
            raise ValueError('For a date range, provide --start, --end (exclusive), and --download')
        first, last = _parse_day(args.start), _parse_day(args.end)
        if first >= last:
            raise ValueError('--end must be after --start')
        import requests
        with requests.Session() as session:
            found = _source_index(session)
            paths, missing_days = [], []
            for offset in range((last - first).days):
                day = first + timedelta(days=offset)
                key = day.strftime('%Y%m%d')
                filename = found.get(key)
                if filename is None:
                    missing_days.append(day.isoformat())
                    print(f'[MISSING FILE] {day} (not a quiet solar-wind day)')
                    continue
                path = _download(session, filename, args.cache)
                paths.append(path)
                print(f'[CDF] {path}', flush=True)
        interval_start, interval_end = first, last
    if not paths:
        raise RuntimeError('No CDF files for the requested range')
    tables = []
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
        tables.append(read_cdf(path, args.variables))
    all_rows = pd.concat(tables, ignore_index=True)
    all_rows['timestamp_utc'] = pd.to_datetime(all_rows['timestamp_utc'], utc=True)
    if interval_start:
        all_rows = all_rows.loc[all_rows.timestamp_utc >= pd.Timestamp(interval_start, tz='UTC')]
    if interval_end:
        all_rows = all_rows.loc[all_rows.timestamp_utc < pd.Timestamp(interval_end, tz='UTC')]
    all_rows = all_rows.sort_values(['variable', 'timestamp_utc']).reset_index(drop=True)
    if all_rows.duplicated(['variable', 'timestamp_utc']).any():
        raise ValueError('Duplicate variable/timestamp across CDF files; refusing ambiguous data')
    if all_rows.empty:
        raise RuntimeError('No measurements found in the requested interval')
    output = args.output or Path('data/isois') / (
        f'summary_rates_{args.start or "local"}_{args.end or "file"}.csv')
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + '.tmp')
    try:
        written = all_rows.copy()
        written['timestamp_utc'] = written['timestamp_utc'].map(_iso)
        written.to_csv(temporary, index=False)
        temporary.replace(output)
    finally:
        temporary.unlink(missing_ok=True)
    manifest = {
        'product': PRODUCT, 'analysis_grade': 'survey_only',
        'not_for_science_analysis': True,
        'source': BASE_URL,
        'time_tags': 'integration midpoints; DELTA_PLUS/MINUS applied',
        'range_start_utc': str(interval_start) if interval_start else None,
        'range_end_exclusive_utc': str(interval_end) if interval_end else None,
        'variables': list(args.variables),
        'source_files': [p.name for p in paths],
        'missing_source_days': missing_days,
        'rows_total': len(all_rows),
        'rows_valid': int(all_rows.valid_measurement.sum()),
        'csv_sha256': _digest(output),
        'warning': 'Summary survey rate product, not calibrated differential flux; do not publish flux correlations from these data.',
    }
    output.with_suffix(output.suffix + '.manifest.json').write_text(
        json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    print(f'[SAVED] {output}: {len(all_rows)} rows, {manifest["rows_valid"]} valid; '
          f'{len(missing_days)} missing source days')
    print('[SURVEY ONLY] NOT calibrated differential flux and NOT a science-grade inference product')


def self_test():
    class FakeCDF:
        def __init__(self, vals, half=None):
            self.vals = vals
            self.half = half if half is not None else [30e9] * 4
        def varget(self, name):
            return {
                'A_H_Rate_TS': self.vals,
                'Epoch': [0, 1, 2, 3],
                'Epoch_DELTA': self.half,
            }[name]
        def varattsget(self, name):
            return {
                'A_H_Rate_TS': {'DEPEND_0': 'Epoch', 'UNITS': 'counts s!E-1!N',
                                'FILLVAL': -1e31, 'VALIDMIN': 0, 'VALIDMAX': 1e6},
                'Epoch': {'DELTA_PLUS_VAR': 'Epoch_DELTA', 'DELTA_MINUS_VAR': 'Epoch_DELTA'},
                'Epoch_DELTA': {'UNITS': 'ns', 'SI_conversion': '1e-9>seconds'},
            }[name]
    def times(_):
        return pd.date_range('2018-11-06', periods=4, freq='min', tz='UTC').tz_localize(None).to_numpy()
    fake = FakeCDF([5.0, -1e31, -2, 1e8])
    frame = extract_rate(fake, 'A_H_Rate_TS', 'psp_isois_l2-summary_20181106_v26.cdf',
                         epoch_converter=times)
    assert list(frame.valid_measurement) == [True, False, False, False]
    assert list(frame.quality_reason)[1:] == ['fill_value', 'below_valid_min', 'above_valid_max']
    assert frame.loc[0, 'integration_start_utc'].startswith('2018-11-05T23:59:30')
    assert frame.loc[0, 'integration_end_utc'].startswith('2018-11-06T00:00:30')
    try:
        extract_rate(FakeCDF(np.zeros((4, 8))), 'A_H_Rate_TS',
                     'psp_isois_l2-summary_20181106_v26.cdf', epoch_converter=times)
    except ValueError as exc:
        assert 'expected 1D' in str(exc)
    else:
        raise AssertionError('Failed to reject unlabelled multidimensional data')
    assert np.isclose(_half_width_seconds([30e9], {'UNITS': 'ns'})[0], 30.0)
    try:
        _half_width_seconds([30], {'UNITS': 'mystery'})
    except ValueError:
        pass
    else:
        raise AssertionError('Unknown integration units were accepted')
    print('PASS: metadata epochs and integration intervals, fill/valid masks, 2D rejection, unknown-unit rejection')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--self-test', action='store_true')
    p.add_argument('--inspect', type=Path, metavar='CDF', help='Print variable/metadata inventory')
    p.add_argument('--file', type=Path, help='Read one locally available summary CDF')
    p.add_argument('--start', help='UTC start date inclusive: YYYY-MM-DD')
    p.add_argument('--end', help='UTC end date exclusive: YYYY-MM-DD')
    p.add_argument('--download', action='store_true', help='Explicitly allow UNH public CDF downloads')
    p.add_argument('--variables', nargs='+', default=list(DEFAULT_VARIABLES),
                   help='Only 1D count-rate variables with per-record integration metadata')
    p.add_argument('--cache', type=Path, default=DEFAULT_CACHE)
    p.add_argument('--output', type=Path)
    a = p.parse_args()
    if a.self_test:
        self_test()
    else:
        try:
            run(a)
        except (ValueError, RuntimeError, OSError, KeyError) as exc:
            p.exit(1, f'[ERROR] {exc}\n')


if __name__ == '__main__':
    main()
