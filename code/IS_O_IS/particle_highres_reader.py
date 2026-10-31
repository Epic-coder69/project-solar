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
PRODUCTS = ('let1-rates10', 'let1-rates60', 'let1-rates300', 'let1-rates3600',
            'let2-rates60', 'het-rates10', 'het-rates60', 'het-rates300',
            'het-rates3600', 'second-rates')
DEFAULT_PRODUCTS = ('let1-rates10', 'let1-rates60', 'het-rates10', 'second-rates')
FILE_RE = re.compile(
    r'^psp_isois-epihi_l2-(' + '|'.join(map(re.escape, PRODUCTS)) +
    r')_(\d{8})_v(\d+)\.cdf$', re.I
)
SCHEMA_VERSION = 2
CHANNEL_GUIDE = {
    'let1-rates10': '10s LET1 H/He energy spectra; He is NOT the original Z6-28 heavy-ion channel',
    'let1-rates60': '60s LET1 species-resolved C, N, O, Ne, Mg, Si, Fe etc.; not identical to the combined Z6-28 summary channel',
    'let1-rates300': '300s LET1 CNO, Ne-Si, Fe-group measurements; not identical to A_Heavy_Rate_TS',
    'let1-rates3600': '3600s LET1 A/B proton and heavy-element spectra; coarse background, not resolved switchbacks',
    'let2-rates60': '60s LET2 C-side spectrum; different detector/view, not a LET1 gap replacement',
    'het-rates10': '10s HET H and electron energy spectra plus He',
    'het-rates60': '60s HET rates/energy spectra where published; coarser event measurements',
    'het-rates300': '300s HET spectra where published; event-scale, not 10s flanks',
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
    url = urljoin(BASE_URL, filename)
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
    match = FILE_RE.fullmatch(path.name)
    if not match or match.group(1).lower() != product:
        raise ValueError(f'Incorrect {product} filename: {path}')
    fingerprint = hashlib.sha256(json.dumps({
        'cdf': _sha256(path), 'catalog': catalog_hash, 'product': product,
        'windows': [[a.isoformat(), b.isoformat()] for a, b in windows]
            if windows is not None else None,
        'schema': SCHEMA_VERSION,
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

    assert _remote_index('psp_isois-epihi_l2-het-rates10_20181106_v2.cdf '
                         'psp_isois-epihi_l2-het-rates10_20181106_v26.cdf',
                         'het-rates10')['20181106'].endswith('_v26.cdf')
    print('PASS: EPI-Hi timing, energy axes, 1s counts, flags, NPZ, network retry and safe Range resume')


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--self-test', action='store_true')
    p.add_argument('--file', type=Path, help='Decode one existing EPI-Hi CDF')
    p.add_argument('--events-csv', type=Path)
    p.add_argument('--products', nargs='+', choices=PRODUCTS, default=list(DEFAULT_PRODUCTS),
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
    if args.self_test:
        self_test()
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
            'schema_version': SCHEMA_VERSION,
            'status': 'complete' if complete else 'partial',
            'fields_catalog': str(args.events_csv) if args.events_csv else None,
            'fields_catalog_sha256': catalog_hash,
            'event_padding_seconds': {'before': before_seconds, 'after': after_seconds},
            'products': ([FILE_RE.fullmatch(args.file.name).group(1)]
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
        match = FILE_RE.fullmatch(args.file.name)
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
                            BASE_URL, timeout=(20, 60),
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
                    m = FILE_RE.fullmatch(cached.name)
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

if __name__ == '__main__':
    main()
