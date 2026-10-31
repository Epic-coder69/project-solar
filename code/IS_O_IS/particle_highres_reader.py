#!/usr/bin/env python3
"""Native-cadence PSP ISʘIS EPI-Hi L2 particle reader.

Keeps detector integration boundaries, energy/sector axes, counts, flux,
uncertainties, viewing geometry and quality flags in separate daily NPZ files.
Default products: LET1 10s (H, He), LET1 60s (heavy elements), HET 10s
(H, He, electrons). The explicitly selected 1s SECOND-RATES product is reduced-bin
quicklook and must not be substituted for full science-grade flux comparisons.
The original hourly summary mode remains unchanged.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import date, timedelta
from pathlib import Path
from urllib.parse import urljoin

import numpy as np
import pandas as pd

from particle_data_reader import (
    _attr, _half_width_seconds, _number, _sha256, _atomic_csv, _atomic_text,
    _utc, validate_fields_events, build_event_windows, merged_intervals,
    required_days,
)

BASE_URL = 'https://spp-isois.sr.unh.edu/data_public/EPIHi/level2/'
PRODUCTS = ('let1-rates10', 'let1-rates60', 'let1-rates300', 'het-rates10', 'second-rates')
DEFAULT_PRODUCTS = ('let1-rates10', 'let1-rates60', 'het-rates10')
FILE_RE = re.compile(
    r'^psp_isois-epihi_l2-(let1-rates10|let1-rates60|let1-rates300|het-rates10|second-rates)_(\d{8})_v(\d+)\.cdf$', re.I
)
SCHEMA_VERSION = 1
CHANNEL_GUIDE = {
    'let1-rates10': '10s LET1 H/He energy spectra; He is NOT the original Z6-28 heavy-ion channel',
    'let1-rates60': '60s LET1 species-resolved C, N, O, Ne, Mg, Si, Fe etc.; not identical to the combined Z6-28 summary channel',
    'let1-rates300': '300s LET1 CNO, Ne-Si, Fe-group measurements; not identical to A_Heavy_Rate_TS',
    'het-rates10': '10s HET H and electron energy spectra plus He',
    'second-rates': '1s reduced-energy-bin H/electron quicklook; not for quantitative science',
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


def decode(cdf, windows=None, epoch_decoder=_cdf_epoch_ns):
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
    details, depends = {}, set()
    for name in names:
        if name == 'Epoch':
            continue
        attrs = cdf.varattsget(name)
        if _attr(attrs, 'DEPEND_0') != 'Epoch':
            continue
        a = _safe_array(cdf.varget(name))
        if a.ndim < 1 or a.shape[0] != len(midpoint):
            raise ValueError(f'{name}: time axis {a.shape} inconsistent with Epoch')
        a = a[keep].copy()
        arrays[name] = a
        # Numeric mask is independent of physical quality flags.
        if np.issubdtype(a.dtype, np.number) and name != 'Quality_Flag' and (
            name.endswith(('_Rate', '_Flux', '_Uncertainty',
                           '_Uncertainty_Plus', '_Uncertainty_Minus'))
            or re.match(r'^(R\d+)?[ABC]_(H|He|Electrons|C|N|O|Fe|Mg|Si)($|_)', name)
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
    if not any(n.endswith('_Flux') for n in arrays):
        candidates = [(n, str(_attr(cdf.varattsget(n), 'DEPEND_0', '')))
                      for n in names if n.endswith(('_Flux', '_Rate'))]
        raise ValueError(
            'No time-aligned flux arrays decoded from CDF. '
            'Inspect source dependencies before assuming data are missing: '
            f'flux/rate candidate names and DEPEND_0={candidates[:24]}; '
            f'Epoch records={len(midpoint)}. '
            'The hourly/10s/60s science products remain available.'
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
        'quality_flags_preserved': 'Quality_Flag' in arrays,
        'quality_flags_interpreted': False,
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
        match = FILE_RE.fullmatch(name)
        if match and match.group(1).lower() == product:
            day, version = match.group(2), int(match.group(3))
            if day not in out or version > out[day][0]:
                out[day] = (version, name)
    return {day: filename for day, (_, filename) in out.items()}


def _download(session, filename, folder):
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / filename
    if path.is_file() and path.stat().st_size:
        return path
    part = folder / (filename + '.part')
    with session.get(urljoin(BASE_URL, filename), stream=True,
                     timeout=(20, 120)) as response:
        response.raise_for_status()
        expected = int(response.headers.get('Content-Length') or 0)
        with part.open('wb') as stream:
            for chunk in response.iter_content(chunk_size=1024*1024):
                if chunk:
                    stream.write(chunk)
    if not part.stat().st_size or (expected and part.stat().st_size != expected):
        raise IOError(f'Incomplete download: {filename}; .part retained')
    import cdflib
    cdf = cdflib.CDF(str(part))
    info = cdf.cdf_info()
    if not (info.zVariables or info.rVariables):
        raise IOError('Downloaded CDF contains no variables')
    part.replace(path)
    return path


def _coverage(windows, starts, ends):
    rows = []
    for w in windows.itertuples(index=False):
        for label, left, right in (
            ('before', w.window_start_utc, w.event_start_utc),
            ('during', w.event_start_utc, w.event_end_exclusive_utc),
            ('after', w.event_end_exclusive_utc, w.window_end_exclusive_utc),
        ):
            a, b = _utc(left).value, _utc(right).value
            if a == b:
                continue
            spans = sorted((max(int(s), a), min(int(e), b))
                           for s, e in zip(starts, ends) if s < b and e > a)
            duration, x0, x1 = 0, None, None
            for s, e in spans:
                if x0 is None:
                    x0, x1 = s, e
                elif s <= x1:
                    x1 = max(x1, e)
                else:
                    duration += x1 - x0
                    x0, x1 = s, e
            if x0 is not None:
                duration += x1 - x0
            rows.append({
                'event_id': w.event_id, 'segment': label,
                'integration_overlap_seconds': duration / 1e9,
                'coverage_fraction_unflagged': min(1.0, duration/(b-a)),
                'overlapping_integrations': len(spans),
            })
    return pd.DataFrame(rows)


def _process(path, product, output_dir, windows, event_df, catalog_hash, resume):
    import cdflib
    match = FILE_RE.fullmatch(path.name)
    if not match or match.group(1).lower() != product:
        raise ValueError(f'Incorrect {product} filename: {path}')
    dest = output_dir / product / (path.stem + '.npz')
    manifest_file = Path(str(dest) + '.manifest.json')
    fingerprint = hashlib.sha256(json.dumps({
        'cdf': _sha256(path), 'catalog': catalog_hash, 'product': product,
        'windows': [[a.isoformat(), b.isoformat()] for a, b in windows]
            if windows is not None else None,
        'schema': SCHEMA_VERSION,
    }, sort_keys=True).encode()).hexdigest()
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
    arrays, info = decode(cdf, windows)
    _save_npz(dest, arrays)
    coverage_path = None
    if event_df is not None:
        coverage_path = dest.with_suffix('.coverage.csv')
        _atomic_csv(_coverage(event_df, arrays['integration_start_utc_ns'],
                              arrays['integration_end_utc_ns']), coverage_path)
    manifest = {
        'schema_version': SCHEMA_VERSION, 'product': product,
        'source_cdf': str(path), 'source_sha256': _sha256(path),
        'output_npz': str(dest), 'output_sha256': _sha256(dest),
        'fingerprint': fingerprint, 'events_sha256': catalog_hash,
        'coverage_csv': str(coverage_path) if coverage_path else None,
        'coverage_sha256': _sha256(coverage_path) if coverage_path else None,
        'timestamp_unit': 'UTC Unix nanoseconds',
        'intervals': 'Epoch midpoint +/- DELTA_MINUS/PLUS',
        'channel_mapping': CHANNEL_GUIDE[product],
        'no_rebinning_or_averaging': True,
        'warning': 'Metadata masks do not replace ISʘIS quality-flag interpretation and count statistics.',
        **info,
    }
    _atomic_text(manifest_file, json.dumps(manifest, indent=2) + '\n')
    return manifest, 'READ'


def self_test():
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
    assert _remote_index('psp_isois-epihi_l2-het-rates10_20181106_v2.cdf '
                         'psp_isois-epihi_l2-het-rates10_20181106_v26.cdf',
                         'het-rates10')['20181106'].endswith('_v26.cdf')
    print('PASS: native EPI-Hi overlap, energy axes, flags, validity, NPZ and version selection')


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--self-test', action='store_true')
    p.add_argument('--file', type=Path, help='Decode one existing EPI-Hi CDF')
    p.add_argument('--events-csv', type=Path)
    p.add_argument('--products', nargs='+', choices=PRODUCTS, default=list(DEFAULT_PRODUCTS),
                   help='Select EPI-Hi products; 1s second-rates requires explicit selection pending metadata validation')
    p.add_argument('--start', help='UTC YYYY-MM-DD inclusive')
    p.add_argument('--end', help='UTC YYYY-MM-DD exclusive')
    p.add_argument('--before-minutes', type=float, default=30)
    p.add_argument('--after-minutes', type=float, default=30)
    p.add_argument('--day-halo', type=int, choices=(0,1), default=1)
    p.add_argument('--cache', type=Path, default=Path('data/isois/epihi_cdf'))
    p.add_argument('--output-dir', type=Path, default=Path('data/isois/highres'))
    p.add_argument('--download', action='store_true')
    p.add_argument('--offline', action='store_true')
    p.add_argument('--resume', action='store_true')
    p.add_argument('--plan-only', action='store_true')
    args = p.parse_args(argv)
    if args.self_test:
        self_test()
        return
    if args.offline and args.download:
        p.error('--offline and --download are mutually exclusive')
    if not np.isfinite(args.before_minutes) or not np.isfinite(args.after_minutes) \
            or min(args.before_minutes, args.after_minutes) < 0:
        p.error('Invalid event padding')
    event_df, windows, catalog_hash = None, None, None
    if args.events_csv:
        if args.start or args.end:
            p.error('--events-csv cannot be combined with --start/--end')
        events = validate_fields_events(args.events_csv)
        event_df = build_event_windows(
            events, before_minutes=args.before_minutes, after_minutes=args.after_minutes
        )
        windows = merged_intervals(event_df)
        catalog_hash = _sha256(args.events_csv)
    paths = []
    missing = []
    if args.file:
        if args.plan_only:
            p.error('--plan-only cannot be used with --file')
        if not args.file.is_file():
            p.error(f'Local file does not exist: {args.file}')
        match = FILE_RE.fullmatch(args.file.name)
        if not match:
            p.error('Unsupported high-resolution EPI-Hi filename')
        paths = [(match.group(1).lower(), args.file)]
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
              f'{len(days)} UTC source days; products: {", ".join(args.products)}')
        if args.plan_only:
            print('[PLAN ONLY] No remote access or output writes')
            return
        import requests
        with requests.Session() as session:
            remote = {}
            if not args.offline:
                r = session.get(BASE_URL, timeout=(15,90))
                r.raise_for_status()
                remote = {name: _remote_index(r.text, name) for name in args.products}
            for product in args.products:
                folder = args.cache / product
                local = {}
                for path in folder.glob('*.cdf'):
                    m = FILE_RE.fullmatch(path.name)
                    if m and m.group(1).lower() == product and path.stat().st_size:
                        day, version = m.group(2), int(m.group(3))
                        if day not in local or version > local[day][0]:
                            local[day] = (version, path)
                for day in days:
                    key = day.strftime('%Y%m%d')
                    if args.offline:
                        path = local.get(key, (None, None))[1]
                    else:
                        name = remote[product].get(key)
                        path = folder/name if name else None
                    if path is None or (not path.exists() and not args.download):
                        missing.append(f'{product} {day}')
                        continue
                    if not path.exists():
                        path = _download(session, path.name, folder)
                    paths.append((product, path))
        if missing:
            print(f'[MISSING] {len(missing)} source product-days: {", ".join(missing[:16])}')
    results = []
    for index, (product, path) in enumerate(paths, 1):
        manifest, action = _process(path, product, args.output_dir, windows,
                                    event_df, catalog_hash, args.resume)
        print(f'[{action} {index}/{len(paths)}] {product}: '
              f'{manifest["retained_time_records"]} integrations; '
              f'{manifest["median_integration_seconds"]}s duration')
        results.append({
            'product': product, 'source_file': str(path), 'npz': manifest['output_npz'],
            'time_records': manifest['retained_time_records'],
        })
    if not results:
        raise RuntimeError('No usable source CDFs: check source-day availability')
    if event_df is not None:
        _atomic_csv(event_df, args.output_dir / 'fields_event_windows.csv')
    _atomic_text(args.output_dir / 'index.json', json.dumps({
        'schema_version': SCHEMA_VERSION,
        'fields_catalog': str(args.events_csv) if args.events_csv else None,
        'fields_catalog_sha256': catalog_hash,
        'files': results,
        'warning': 'No temporal averaging. Quality and particle counting statistics still need analysis.',
    }, indent=2) + '\n')
    print(f'[SAVED] {len(results)} native-cadence partitions; index: {args.output_dir / "index.json"}')


if __name__ == '__main__':
    main()
