#!/usr/bin/env python3
"""Project Solar: full-cadence magnetic-event analysis plus legacy benchmark.

Laker et al. (2022) event thresholds (45/30 deg, 10 s, merge <20 s).
References:
https://academic.oup.com/mnras/article/517/1/1001/6692877
SWEAP SPC Level 3: PSP_SWP_SPC_L3I, vp1_fit_RTN, sc_pos_HCI, general_flag.

--reference local: 3-hour component-wise median *direction*, but independent
                   rolling median of scalar |B| for magnitude deviation.
--reference parker-spc: APPROXIMATE Parker spiral using SPC radial proton speed,
                   HCI position, assumed source surface, and local-sector polarity.
                   This is not a reproduction of the published SPAN-i + SPC
                   fallback and filtered velocity method.

Default --mode full consumes the reader's --windows-csv output and analyzes
original full-cadence FIELDS data. Detection boundaries and plot traces are
at one-second resolution; native samples are revisited for physical metrics.
--mode legacy preserves the older 1min/4cyc/SPC benchmark workflow.
Never labels any magnetic deflection as a confirmed physical switchback.
"""
import argparse
from datetime import timedelta
from functools import lru_cache
from pathlib import Path
from urllib.parse import urlparse
import re
import numpy as np
import pandas as pd
import requests
import cdflib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from magnetic_data_reader import get_day, read_high, full_cadence_seconds, read_full_interval

DATA = Path('data')
FIG = Path('figures')
SPC_CACHE = DATA / 'sweap_spc'
SPC_API = ('https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/sp_phys/'
           'datasets/PSP_SWP_SPC_L3I/orig_data')
TRIGGER_DEG, EDGE_DEG = 45., 30.
MIN_SECONDS, MERGE_SECONDS = 10, 20
SOLAR_RADIUS_KM = 695700.
OMEGA_SUN_RAD_S = 2 * np.pi / (25.38 * 86400.)

# Frozen November 2018 exploratory intervals from the original survey.
# These are candidate-selection windows, NOT confirmed physical switchbacks.
# The original saved start/end values were reconstructed during cleanup;
# check against the user's original local catalog before publication-grade reruns.
FROZEN_SURVEY_WINDOWS = (
    ("2018-11-27T13:22:30", "2018-11-27T13:31:30", "strong"),
    ("2018-11-20T21:58:30", "2018-11-20T22:19:30", "strong"),
    ("2018-11-07T03:32:30", "2018-11-07T03:54:30", "strong"),
    ("2018-11-01T18:42:30", "2018-11-01T18:43:30", "strong"),
    ("2018-11-27T19:38:30", "2018-11-27T19:41:30", "strong"),
    ("2018-11-08T09:44:30", "2018-11-08T09:49:30", "strong"),
    ("2018-11-27T04:33:30", "2018-11-27T04:46:30", "strong"),
    ("2018-11-05T13:50:30", "2018-11-05T13:55:30", "strong"),
    ("2018-11-25T21:29:30", "2018-11-25T21:51:30", "strong"),
    ("2018-11-29T07:47:30", "2018-11-29T07:57:30", "strong"),
    ("2018-11-08T09:38:30", "2018-11-08T09:41:30", "ambiguous"),
    ("2018-11-26T20:25:30", "2018-11-26T20:26:30", "ambiguous"),
    ("2018-11-28T07:17:30", "2018-11-28T07:21:30", "ambiguous"),
    ("2018-11-27T14:50:30", "2018-11-27T14:54:30", "ambiguous"),
    ("2018-11-20T08:51:30", "2018-11-20T08:52:30", "ambiguous"),
    ("2018-11-21T05:50:30", "2018-11-21T05:50:30", "control"),
    ("2018-11-12T04:20:30", "2018-11-12T04:20:30", "control"),
    ("2018-11-06T12:50:30", "2018-11-06T12:50:30", "control"),
    ("2018-11-26T22:20:30", "2018-11-26T22:20:30", "control"),
    ("2018-11-02T02:05:30", "2018-11-02T02:05:30", "control"),
)


def default_windows():
    """Return the 20 frozen survey intervals without requiring a CSV file."""
    return pd.DataFrame(FROZEN_SURVEY_WINDOWS, columns=["start", "end", "category"])

@lru_cache(maxsize=5)
def high_for_day(day):
    return read_high(day)

@lru_cache(maxsize=8)
def one_min_day(day):
    return get_day(day)


def date_sequence(left, right, pad_days=0):
    return [d.date() for d in pd.date_range(
        left.normalize()-pd.Timedelta(days=pad_days),
        right.normalize()+pd.Timedelta(days=pad_days), freq='D')]


def scalar_and_direction_reference(left, right):
    chunks = []
    for d in date_sequence(left, right, pad_days=1):
        try:
            chunks.append(one_min_day(d)[['BR','BT','BN','Bmag']])
        except Exception as exc:
            print(f'[FIELDS one-minute unavailable] {d}: {exc}')
    if not chunks:
        raise RuntimeError('No one-minute measurements for background')
    raw = pd.concat(chunks).sort_index()
    raw = raw[~raw.index.duplicated(keep='last')]
    # A local directional baseline is a vector of independently smoothed components.
    direction = raw[['BR','BT','BN']].rolling('181min', center=True, min_periods=60).median()
    # CRITICAL: scalar magnitude baseline must smooth original scalar |B|, NOT |direction|.
    magnitude = raw['Bmag'].rolling('181min', center=True, min_periods=60).median()
    return direction, magnitude


def latest_version(name):
    hit = re.search(r'_v(\d+)\.cdf$', name, re.I)
    return int(hit.group(1)) if hit else -1


def spc_file(day):
    SPC_CACHE.mkdir(parents=True, exist_ok=True)
    marker = f'{day:%Y%m%d}'
    prior = [p for p in SPC_CACHE.glob(f'*{marker}*_v*.cdf') if 'spc_l3i' in p.name.lower()]
    if prior:
        return max(prior, key=lambda p: latest_version(p.name))
    stop = day + timedelta(days=1)
    url = f'{SPC_API}/{day:%Y%m%d}T000000Z,{stop:%Y%m%d}T000000Z/'
    resp = requests.get(url, headers={'Accept':'application/json'}, timeout=90)
    resp.raise_for_status()
    descriptions = resp.json().get('FileDescription', [])
    if isinstance(descriptions, dict):
        descriptions = [descriptions]
    files = [(Path(urlparse(item['Name']).path).name, item['Name']) for item in descriptions
             if item.get('Name', '').lower().endswith('.cdf') and marker in item['Name']]
    if not files:
        raise RuntimeError(f'No SWEAP SPC L3I file for {day}')
    name, remote = max(files, key=lambda p: latest_version(p[0]))
    target = SPC_CACHE / name
    partial = SPC_CACHE / (name + '.part')
    try:
        with requests.get(remote, stream=True, timeout=(30, 180)) as r:
            r.raise_for_status()
            with partial.open('wb') as output:
                for chunk in r.iter_content(1024*1024):
                    if chunk:
                        output.write(chunk)
        if partial.stat().st_size < 1024:
            raise ValueError('SPC download unexpectedly small')
        partial.replace(target)
    finally:
        partial.unlink(missing_ok=True)
    print('[SPC downloaded]', target)
    return target


@lru_cache(maxsize=5)
def spc_day(day):
    cdf = cdflib.CDF(str(spc_file(day)))
    names = {n.lower(): n for n in list(cdf.cdf_info().zVariables) + list(cdf.cdf_info().rVariables)}
    needed = ('vp1_fit_rtn', 'sc_pos_hci', 'general_flag')
    missing = [name for name in needed if name not in names]
    if missing:
        raise RuntimeError(f'SPC variables missing: {missing}; available: {list(names)}')
    vector_var = names['vp1_fit_rtn']
    epoch = cdf.varattsget(vector_var).get('DEPEND_0')
    if isinstance(epoch, (list, tuple, np.ndarray)):
        epoch = epoch[0]
    if isinstance(epoch, bytes):
        epoch = epoch.decode()
    if not epoch:
        epoch = names.get('epoch')
    if epoch is None:
        raise RuntimeError('SPC epoch variable not found')
    idx = pd.to_datetime(cdflib.cdfepoch.to_datetime(cdf.varget(epoch)))
    vel = np.asarray(cdf.varget(vector_var), dtype=float)
    pos = np.asarray(cdf.varget(names['sc_pos_hci']), dtype=float)
    quality = np.asarray(cdf.varget(names['general_flag'])).ravel()
    if len(idx) != len(vel) or vel.ndim != 2 or vel.shape[1] != 3 or pos.shape != vel.shape or len(quality) != len(idx):
        raise RuntimeError('Unexpected shape for SWEAP velocity, position, or quality flags')
    frame = pd.DataFrame({'vr':vel[:,0], 'x':pos[:,0], 'y':pos[:,1], 'z':pos[:,2], 'flag':quality}, index=idx)
    # Only measurements with NASA nominal (zero) quality flag are used.
    frame.loc[frame.flag != 0, ['vr','x','y','z']] = np.nan
    for col in ('vr','x','y','z'):
        frame.loc[(~np.isfinite(frame[col])) | (frame[col].abs() >= 1e20), col] = np.nan
    frame.loc[~frame.vr.between(100, 1500), 'vr'] = np.nan
    frame['r'] = np.sqrt(frame.x**2 + frame.y**2 + frame.z**2)
    frame['sin_colat'] = np.sqrt(frame.x**2 + frame.y**2) / frame.r
    frame = frame[~frame.index.duplicated(keep='last')].sort_index()
    # Reduce high-cadence SPC moments before computing low-frequency reference.
    return frame[['vr','r','sin_colat']].resample('1min').median()


def parker_reference(left, right, direction):
    chunks = []
    for d in date_sequence(left, right, pad_days=1):
        try:
            chunks.append(spc_day(d))
        except Exception as exc:
            print(f'[SPC unavailable] {d}: {exc}')
    if not chunks:
        raise RuntimeError('No SPC data available; use --reference local or inspect SPC archive')
    plasma = pd.concat(chunks).sort_index()
    plasma = plasma[~plasma.index.duplicated(keep='last')]
    # Exploratory smoothing. Laker et al. used a 2nd-order Butterworth 3-hour filter
    # and generally SPAN-i core fits with SPC as gap filler. This is NOT identical.
    params = plasma.rolling('181min', center=True, min_periods=60).median()
    ratio = -OMEGA_SUN_RAD_S * (params.r - 2.5*SOLAR_RADIUS_KM) * params.sin_colat / params.vr
    ratio[(params.r <= 2.5*SOLAR_RADIUS_KM) | (~np.isfinite(ratio))] = np.nan
    # Use stable local magnetic sector polarity only; solar-wind speed determines winding.
    polarity = np.sign(direction.BR).replace(0, np.nan)
    ratio = ratio.reindex(direction.index, method='nearest', tolerance=pd.Timedelta(seconds=35))
    ref = pd.DataFrame({'BR':polarity, 'BT':polarity*ratio, 'BN':0.0}, index=direction.index)
    ref['BN'] = ref['BN'].where(ref.BR.notna() & ref.BT.notna(), np.nan)
    return ref


def get_window(seed, reference='local', margin_minutes=40):
    center = pd.Timestamp(seed['start'])
    stop = pd.Timestamp(seed['end']) if pd.notna(seed.get('end')) else center
    if stop < center:
        stop = center
    left = center - pd.Timedelta(minutes=margin_minutes)
    right = stop + pd.Timedelta(minutes=margin_minutes)
    chunks = []
    for d in date_sequence(left, right):
        try:
            chunks.append(high_for_day(d))
        except Exception as exc:
            print(f'[HIGH RES unavailable] {d}: {exc}')
    if not chunks:
        raise RuntimeError(f'No FIELDS 4cyc data for {center}')
    high = pd.concat(chunks).sort_index()
    high = high[~high.index.duplicated(keep='last')]
    high = high.loc[left:right, ['BR','BT','BN']]
    if high.empty or high.notna().all(axis=1).sum() < 20:
        raise RuntimeError('Insufficient usable high-resolution measurements')
    b = high.resample('1s').median()
    b['Bmag'] = np.sqrt((b[['BR','BT','BN']]**2).sum(axis=1, min_count=3))
    local_direction, local_magnitude = scalar_and_direction_reference(left, right)
    ref_minute = local_direction if reference == 'local' else parker_reference(left, right, local_direction)
    ref = ref_minute.reindex(b.index, method='nearest', tolerance=pd.Timedelta(seconds=35))
    magnitude = local_magnitude.reindex(b.index, method='nearest', tolerance=pd.Timedelta(seconds=35))
    matched = (b[['BR','BT','BN']].notna().all(axis=1) & ref[['BR','BT','BN']].notna().all(axis=1))
    if int(matched.sum()) < 20 or matched.mean() < .1:
        raise RuntimeError(f'Insufficient matched {reference} reference coverage ({matched.mean():.1%})')
    return b, ref, magnitude, left, right


def deflection(b, ref):
    a = b[['BR','BT','BN']].to_numpy(float)
    r = ref[['BR','BT','BN']].to_numpy(float)
    norms = np.linalg.norm(a,axis=1)*np.linalg.norm(r,axis=1)
    valid = np.isfinite(norms)&(norms>1e-8)&np.isfinite(a).all(axis=1)&np.isfinite(r).all(axis=1)
    cos = np.full(len(b),np.nan)
    cos[valid] = np.sum(a[valid]*r[valid],axis=1)/norms[valid]
    return pd.Series(np.degrees(np.arccos(np.clip(cos,-1,1))), index=b.index)


def find_events(angle):
    """Angle threshold detection; a NaN splits rather than bridging a data gap."""
    hit = angle.ge(EDGE_DEG).fillna(False).to_numpy(bool)
    valid = angle.notna().to_numpy(bool)
    locations = np.flatnonzero(hit)
    if not len(locations):
        return []
    groups = np.split(locations,np.where(np.diff(locations)>1)[0]+1)
    spans = [(int(g[0]),int(g[-1])) for g in groups]
    merged = [spans[0]]
    for start,end in spans[1:]:
        before, old_end = merged[-1]
        gap = (angle.index[start]-angle.index[old_end]).total_seconds()
        if gap < MERGE_SECONDS and valid[old_end:start+1].all():
            merged[-1] = before,end
        else:
            merged.append((start,end))
    result = []
    for i,j in merged:
        observed = (angle.index[j]-angle.index[i]).total_seconds()+1
        peak = float(angle.iloc[i:j+1].max())
        if observed < MIN_SECONDS or peak < TRIGGER_DEG:
            continue
        edge_start = i==0 or (not valid[i-1]) or (angle.iloc[i-1]>=EDGE_DEG)
        edge_end = j==len(angle)-1 or (not valid[j+1]) or (angle.iloc[j+1]>=EDGE_DEG)
        result.append((i,j,observed,peak,edge_start,edge_end))
    return result


def evaluate_seed(seed, number, reference='local', plot=True):
    b, ref, mag_ref, requested_start, requested_end = get_window(seed,reference)
    angle = deflection(b,ref)
    events = []
    center = pd.Timestamp(seed['start'])
    stop = pd.Timestamp(seed['end']) if pd.notna(seed.get('end')) else center
    cat = str(seed.get('category','unknown'))
    for i,j,observed,peak,clip_left,clip_right in find_events(angle):
        start,end = b.index[i],b.index[j]
        if end < center - pd.Timedelta(minutes=15) or start > stop + pd.Timedelta(minutes=15):
            continue
        section = b.iloc[i:j+1]
        reference_slice = ref.iloc[i:j+1]
        ratio = mag_ref.iloc[i:j+1]
        valid = (ratio > 0.01) & section.Bmag.notna()
        deviations = ((section.Bmag[valid] / ratio[valid])-1).abs()*100
        reversal = ((section.BR * reference_slice.BR) < 0).any()
        censored = bool(clip_left or clip_right)
        if not len(deviations):
            median_delta = max_delta = np.nan
        else:
            median_delta, max_delta = float(deviations.median()), float(deviations.max())
        events.append(dict(
            seed_utc=center.isoformat(), seed_category=cat,
            event_start_utc=start.isoformat(), event_end_utc=end.isoformat(),
            observed_span_seconds=observed,
            duration_seconds=(np.nan if censored else observed),
            left_censored=bool(clip_left), right_censored=bool(clip_right),
            maximum_deflection_deg=round(peak,1),
            median_abs_Bmag_deviation_pct=round(median_delta,2),
            max_abs_Bmag_deviation_pct=round(max_delta,2),
            radial_reversal=bool(reversal),
            reference=reference, status='unvalidated_magnetic_rotation',
        ))
    image_path = ''
    if plot:
        FIG.mkdir(exist_ok=True)
        image_path = str(FIG/f'corrected_{reference}_{cat}_{center:%Y%m%d_%H%M}_{number}.png')
        fig, ax = plt.subplots(3,1,figsize=(12,9),sharex=True,constrained_layout=True)
        for col in ('BR','BT','BN'):
            ax[0].plot(b.index,b[col],linewidth=.7,label=col)
        ax[0].legend(ncol=3); ax[0].set_ylabel('RTN field (nT)')
        ax[1].plot(b.index,b.Bmag,color='black',linewidth=.8,label='Measured |B|')
        ax[1].plot(b.index,mag_ref,color='tab:purple',lw=.9,ls='--',label='Scalar |B| background')
        ax[1].legend(); ax[1].set_ylabel('|B| (nT)')
        ax[2].plot(angle.index,angle,lw=.8)
        ax[2].axhline(EDGE_DEG,color='gray',ls='--',label='30° boundary')
        ax[2].axhline(TRIGGER_DEG,color='red',ls='--',label='45° trigger')
        ax[2].legend(ncol=2); ax[2].set_ylabel('Deflection (degrees)'); ax[2].set_xlabel('UTC')
        for event in events:
            l = pd.Timestamp(event['event_start_utc']); r = pd.Timestamp(event['event_end_utc'])
            for axis in ax:
                axis.axvspan(l,r,color='orange',alpha=.15,
                            hatch='///' if (event['left_censored'] or event['right_censored']) else None)
        for axis in ax:
            axis.grid(alpha=.2)
        fig.suptitle(f'Corrected benchmark ({reference}) — {cat} — {center:%Y-%m-%d %H:%M} UTC')
        fig.savefig(image_path,dpi=150); plt.close(fig)
    return events,image_path


def deduplicate(events):
    """Exact matches and substantial overlaps from different survey seed windows."""
    if not events:
        return []
    x = pd.DataFrame(events).copy()
    x['_start'] = pd.to_datetime(x.event_start_utc)
    x['_end'] = pd.to_datetime(x.event_end_utc)
    x = x.sort_values(['_start','_end'])
    groups = []
    for row in x.to_dict('records'):
        found = None
        for group in groups:
            first = group[0]
            # Events found in *different* windows; at least half of the
            # shorter event must overlap. Overlapping physics can remain distinct.
            overlap = (min(row['_end'],first['_end']) - max(row['_start'],first['_start'])).total_seconds()+1
            shorter = min((row['_end']-row['_start']).total_seconds()+1,
                          (first['_end']-first['_start']).total_seconds()+1)
            if row['reference']==first['reference'] and shorter>0 and overlap/shorter >= .5:
                found = group
                break
        if found is None:
            groups.append([row])
        else:
            found.append(row)
    unique=[]
    for group in groups:
        # A complete event is preferable to censored measurements of same interval.
        group.sort(key=lambda r:(bool(r['left_censored'] or r['right_censored']),
                                 -(r['_end']-r['_start']).total_seconds()))
        keep=group[0].copy()
        keep.pop('_start',None);keep.pop('_end',None)
        keep['duplicate_source_count']=len(group)
        keep['source_seed_utcs']=';'.join(sorted(set(r['seed_utc'] for r in group)))
        unique.append(keep)
    return sorted(unique,key=lambda r:r['event_start_utc'])


def self_test():
    t=pd.date_range('2018-11-06',periods=200,freq='s')
    a=pd.Series(np.zeros(len(t)),index=t);a.iloc[40:80]=120
    ev=find_events(a)
    assert len(ev)==1 and ev[0][2]==40 and not ev[0][4] and not ev[0][5]
    a.iloc[:]=0;a.iloc[40:45]=120
    assert len(find_events(a))==0
    a.iloc[:]=0;a.iloc[40:55]=120;a.iloc[65:78]=120
    assert len(find_events(a))==1
    a.iloc[:]=0;a.iloc[:30]=120
    assert find_events(a)[0][4] is True
    a.iloc[:]=0;a.iloc[170:]=120
    assert find_events(a)[0][5] is True
    a.iloc[:]=0;a.iloc[40:55]=120;a.iloc[56:66]=120;a.iloc[55]=np.nan
    assert len(find_events(a))==2, 'Must not merge through missing sample'
    # Scalar versus vector median differ for rotating fields.
    vectors=pd.DataFrame({'BR':[10,-10,10,-10], 'BT':[0,0,0,0], 'BN':[0,0,0,0]})
    scalar=np.sqrt((vectors**2).sum(axis=1))
    assert scalar.median()==10
    assert (vectors.median()**2).sum()**.5==0 # demonstrates original denominator failure
    base={'seed_utc':'2018-11-06T00:00:00','seed_category':'strong',
          'event_start_utc':'2018-11-06T00:00:10','event_end_utc':'2018-11-06T00:00:20',
          'left_censored':False,'right_censored':False,'reference':'local'}
    another={**base,'seed_utc':'2018-11-06T00:01:00'}
    assert len(deduplicate([base,another]))==1
    print('PASS: thresholds, duration, NaN gaps, boundary flags, and dedup')



# ---------------------------------------------------------------------------
# Full-cadence candidate analysis (new path). The legacy 4cyc benchmark above
# remains available with --mode legacy for reproducibility.
# ---------------------------------------------------------------------------

FULL_EVENT_COLUMNS = [
    'seed_utc', 'seed_category', 'event_start_utc', 'event_end_utc',
    'observed_span_seconds', 'duration_seconds', 'left_censored',
    'right_censored', 'maximum_deflection_deg',
    'native_peak_deflection_deg', 'native_valid_samples',
    'native_total_samples', 'median_abs_Bmag_deviation_pct_approx',
    'max_abs_Bmag_deviation_pct', 'radial_reversal',
    'reference', 'status', 'duplicate_source_count', 'source_seed_utcs'
]


@lru_cache(maxsize=4)
def _full_seconds_for_day(day):
    """Screening-rate summary derived ONLY from original full-cadence CDFs.

    Caching one-second daily vectors avoids reading a whole day of NASA data
    again for every candidate in the same daily scan.
    """
    return full_cadence_seconds(day)


def _full_window(seed):
    """One candidate's data, with 100-minute background on both sides.

    The candidate's start/end are observation context, not event boundaries.
    Missing-day data are not silently interpolated.
    """
    left, right = pd.Timestamp(seed['start']), pd.Timestamp(seed['end'])
    if left.tzinfo is not None:
        left = left.tz_convert('UTC').tz_localize(None)
    if right.tzinfo is not None:
        right = right.tz_convert('UTC').tz_localize(None)
    if right <= left:
        raise ValueError('Candidate end must follow start')
    halo = pd.Timedelta(minutes=100)
    background_left, background_right = left-halo, right+halo
    pieces = []
    unavailable = []
    for day in pd.date_range(background_left.normalize(),
                             background_right.normalize(), freq='D'):
        try:
            piece = _full_seconds_for_day(day.date())
        except (FileNotFoundError, RuntimeError) as exc:
            unavailable.append(f'{day.date()}: {exc}')
            continue
        pieces.append(piece.loc[(piece.index >= background_left) &
                                (piece.index <= background_right)])
    if not pieces:
        raise RuntimeError('No full-cadence FIELDS measurements: ' +
                           '; '.join(unavailable))
    seconds = pd.concat(pieces).sort_index()
    seconds = seconds[~seconds.index.duplicated(keep='last')]
    # Reinsert missing one-second bins so find_events cannot bridge a data gap.
    seconds = seconds.resample('1s').median()
    sample = seconds.loc[(seconds.index >= left) &
                         (seconds.index <= right)].copy()
    if sample.empty or sample[['BR', 'BT', 'BN']].notna().all(axis=1).sum() < 10:
        raise RuntimeError('Too few valid one-second magnetic vectors in candidate')

    minute = seconds[['BR', 'BT', 'BN']].resample('1min').median()
    direction = minute.rolling('181min', center=True, min_periods=30).median()
    # Scalar baseline is independent of smoothed vector length.
    scalar = seconds['Bmag'].resample('1min').median().rolling(
        '181min', center=True, min_periods=30
    ).median()
    reference = direction.reindex(
        sample.index, method='nearest', tolerance=pd.Timedelta(seconds=35))
    magnitude = scalar.reindex(
        sample.index, method='nearest', tolerance=pd.Timedelta(seconds=35))
    sample['Bmag'] = np.sqrt(
        (sample[['BR', 'BT', 'BN']]**2).sum(axis=1, min_count=3))
    return sample, reference, magnitude, left, right, unavailable


def _native_event_statistics(left, right, reference, mag_reference,
                             native_stream=None):
    """Revisit ORIGINAL magnetic records; never keep a full event in RAM.

    Event boundaries are 1-second-resolution threshold crossings. Peak angle,
    radial reversal and magnitude deviations use native-cadence observations.
    The median deviation is quantized to 0.1 percentage points (streaming
    histogram), hence the _approx name; its upper bin is >=500%.
    """
    if native_stream is None:
        native_stream = read_full_interval(left, right + pd.Timedelta(seconds=1))
    bins = np.zeros(5001, dtype=np.int64)
    total, matched, magnitude_count = 0, 0, 0
    peak, max_delta, reversal = np.nan, np.nan, False
    stop = right + pd.Timedelta(seconds=1)
    for chunk in native_stream:
        chunk = chunk.loc[(chunk.index >= left) & (chunk.index < stop)]
        if chunk.empty:
            continue
        total += len(chunk)
        seconds = chunk.index.floor('s')
        bg = reference.reindex(seconds)[['BR','BT','BN']].to_numpy(dtype=float)
        mag_bg = mag_reference.reindex(seconds).to_numpy(dtype=float)
        field = chunk[['BR','BT','BN']].to_numpy(dtype=float)
        denom = np.linalg.norm(bg, axis=1) * np.linalg.norm(field, axis=1)
        valid = (np.isfinite(bg).all(axis=1) &
                 np.isfinite(field).all(axis=1) &
                 np.isfinite(denom) & (denom > 1e-8))
        matched += int(valid.sum())
        if valid.any():
            cos = np.sum(field[valid]*bg[valid], axis=1)/denom[valid]
            angles = np.degrees(np.arccos(np.clip(cos, -1, 1)))
            value = float(angles.max())
            peak = value if not np.isfinite(peak) else max(peak, value)
            reversal |= bool(np.any(field[valid, 0]*bg[valid, 0] < 0))
        scalar = chunk['Bmag'].to_numpy(dtype=float)
        good = (valid & np.isfinite(scalar) & np.isfinite(mag_bg) &
                (mag_bg > 0.01))
        if good.any():
            deviation = np.abs(scalar[good]/mag_bg[good]-1)*100
            magnitude_count += len(deviation)
            top = float(deviation.max())
            max_delta = top if not np.isfinite(max_delta) else max(max_delta, top)
            positions = np.minimum((deviation*10).astype(int), 5000)
            bins += np.bincount(positions, minlength=5001)
    if magnitude_count:
        median_bin = int(np.searchsorted(np.cumsum(bins),
                                         (magnitude_count+1)//2))
        median_delta = (median_bin + 0.5)/10
    else:
        median_delta = np.nan
    return dict(
        native_total_samples=total,
        native_valid_samples=matched,
        native_peak_deflection_deg=(round(peak, 1)
                                    if np.isfinite(peak) else np.nan),
        median_abs_Bmag_deviation_pct_approx=(
            round(median_delta, 2) if np.isfinite(median_delta) else np.nan),
        max_abs_Bmag_deviation_pct=(
            round(max_delta, 2) if np.isfinite(max_delta) else np.nan),
        radial_reversal=bool(reversal)
    )


def _seed_core(seed, left, right):
    """Core limits choose RELEVANT detections; they never define an event."""
    a = pd.Timestamp(seed.get('core_start', left))
    b = pd.Timestamp(seed.get('core_end', right))
    if a.tzinfo is not None:
        a = a.tz_convert('UTC').tz_localize(None)
    if b.tzinfo is not None:
        b = b.tz_convert('UTC').tz_localize(None)
    if b < a:
        raise ValueError('core_end precedes core_start')
    return a, b


def _plot_full_candidate(seed, number, sample, reference, magnitude,
                         angle, events):
    FIG.mkdir(parents=True, exist_ok=True)
    path = FIG / f'full_analysis_{pd.Timestamp(seed["start"]):%Y%m%d_%H%M%S}_{number}.png'
    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True,
                             constrained_layout=True)
    for name in ('BR', 'BT', 'BN'):
        axes[0].plot(sample.index, sample[name], lw=.7, label=name)
    axes[0].set_ylabel('RTN magnetic field (nT)')
    axes[0].legend(ncol=3)
    axes[1].plot(sample.index, sample['Bmag'], lw=.8, label='Measured |B| (1s)')
    axes[1].plot(magnitude.index, magnitude, lw=.8, ls='--',
                 label='Local scalar reference')
    axes[1].set_ylabel('|B| (nT)')
    axes[1].legend()
    axes[2].plot(angle.index, angle, lw=.8)
    axes[2].axhline(EDGE_DEG, color='grey', ls='--', label='30° boundary')
    axes[2].axhline(TRIGGER_DEG, color='red', ls='--', label='45° trigger')
    axes[2].set_ylabel('Deflection (degrees)')
    axes[2].set_xlabel('UTC')
    axes[2].legend()
    for e in events:
        for ax in axes:
            ax.axvspan(pd.Timestamp(e['event_start_utc']),
                       pd.Timestamp(e['event_end_utc']),
                       color='orange', alpha=.15)
    for ax in axes:
        ax.grid(alpha=.2)
    fig.suptitle('FIELDS full-cadence source, one-second event boundaries '
                 '(unvalidated)')
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return str(path)


def evaluate_full_seed(seed, number, plot=True):
    """Analyze candidate with local reference and native-cadence measurements."""
    sample, reference, mag_ref, left, right, unavailable = _full_window(seed)
    angle = deflection(sample, reference)
    core_left, core_right = _seed_core(seed, left, right)
    events = []
    for i, j, duration, angle_peak, clip_left, clip_right in find_events(angle):
        event_left, event_right = sample.index[i], sample.index[j]
        if event_right < core_left or event_left > core_right:
            continue
        native = _native_event_statistics(
            event_left, event_right, reference, mag_ref)
        censored = bool(clip_left or clip_right)
        events.append(dict(
            seed_utc=left.isoformat(),
            seed_category=str(seed.get('category', 'unknown')),
            event_start_utc=event_left.isoformat(),
            event_end_utc=event_right.isoformat(),
            observed_span_seconds=duration,
            duration_seconds=(np.nan if censored else duration),
            left_censored=bool(clip_left),
            right_censored=bool(clip_right),
            maximum_deflection_deg=round(angle_peak, 1),
            **native, reference='full_local_181min',
            status='unvalidated_magnetic_rotation'
        ))
    figure = (_plot_full_candidate(seed, number, sample, reference, mag_ref,
                                   angle, events) if plot else '')
    return events, figure, dict(
        valid_seconds=int(sample[['BR','BT','BN']].notna().all(axis=1).sum()),
        total_seconds=len(sample),
        matched_reference_seconds=int(angle.notna().sum()),
        missing_background_days='; '.join(unavailable)
    )


def full_self_test():
    """No NASA download: test native vector metrics and core filtering."""
    idx = pd.date_range('2018-11-02 13:00:00', periods=180, freq='s')
    measured = pd.DataFrame({
        'BR': np.where((np.arange(180)>=65)&(np.arange(180)<90), -80., 80.),
        'BT': 0., 'BN': 0.
    }, index=idx)
    measured['Bmag'] = 80.
    ref = pd.DataFrame({'BR':80., 'BT':0., 'BN':0.}, index=idx)
    mag = pd.Series(80., index=idx)
    ang = deflection(measured, ref)
    found = find_events(ang)
    assert len(found) == 1 and found[0][2] == 25
    raw_idx = pd.date_range(idx[65], periods=75, freq='333ms')
    original = pd.DataFrame({'BR':-80.,'BT':0.,'BN':0.,'Bmag':80.},
                            index=raw_idx)
    stats = _native_event_statistics(idx[65], idx[89], ref, mag,
                                     native_stream=iter([original]))
    assert stats['native_valid_samples'] == 75
    assert stats['native_peak_deflection_deg'] == 180.0
    assert stats['radial_reversal'] is True
    assert stats['median_abs_Bmag_deviation_pct_approx'] < .1
    a, b = _seed_core({'core_start':idx[65], 'core_end':idx[89]},
                      idx[0], idx[-1])
    assert (a, b) == (idx[65], idx[89])
    # Exercise the new candidate -> 1-second background path without NASA.
    from unittest.mock import patch
    synthetic_idx = pd.date_range('2018-11-02 10:00:00',
                                  periods=21600, freq='s')
    syn = pd.DataFrame({'BR':80.,'BT':0.,'BN':0.,'Bmag':80.},
                       index=synthetic_idx)
    syn.loc[idx[65]:idx[89], 'BR'] = -80.
    with patch(__name__ + '._full_seconds_for_day', return_value=syn):
        sliced, baseline, scalar, *_ = _full_window({
            'start':'2018-11-02 13:00:00',
            'end':'2018-11-02 13:03:00',
            'category':'screened_full',
            'core_start':idx[65], 'core_end':idx[89]
        })
    assert len(find_events(deflection(sliced, baseline))) == 1
    assert float(scalar.loc[idx[65]]) == 80.
    print('PASS: native-cadence physical metrics, 1s event segmentation, '
          'and candidate-core handling')


def run_full_analysis(args, parser):
    """New CSV-to-physical-analysis pipeline; never overwrites legacy outputs."""
    if args.windows_csv is None:
        parser.error('Full-cadence mode needs --windows-csv from the reader')
    if args.reference != 'local':
        parser.error('Full-cadence mode supports --reference local only; '
                     'use --mode legacy for the older Parker-SPC benchmark')
    if not args.windows_csv.is_file():
        parser.error(f'Candidate CSV not found: {args.windows_csv}')
    seeds = pd.read_csv(args.windows_csv)
    needed = {'start', 'end', 'category'}
    if not needed.issubset(seeds.columns):
        parser.error('CSV requires columns start,end,category')
    if seeds.empty:
        print('No candidate rows in the CSV; no physical events to analyze.')
        return
    seeds = seeds.sort_values('start')
    if args.top:
        seeds = seeds.head(args.top)
    events, audit = [], []
    DATA.mkdir(parents=True, exist_ok=True)
    # Output names follow the input catalog, so two different observation
    # days do not overwrite each other's research results.
    catalog_name = re.sub(r'[^A-Za-z0-9_-]+', '_', args.windows_csv.stem)
    suffix = (f'_top{args.top}' if args.top else '')
    prefix = f'full_analysis_{catalog_name}{suffix}'
    output = DATA / f'{prefix}_events.csv'
    audit_output = DATA / f'{prefix}_audit.csv'

    def save_progress():
        # Tiny catalogs are atomically replaced after each candidate; safe to
        # open from another terminal while a long analysis continues.
        unique_now = deduplicate(events)
        for path, frame in (
            (output, pd.DataFrame(unique_now, columns=FULL_EVENT_COLUMNS)),
            (audit_output, pd.DataFrame(audit)),
        ):
            temporary = path.with_name(path.name + '.tmp')
            try:
                frame.to_csv(temporary, index=False)
                temporary.replace(path)
            finally:
                temporary.unlink(missing_ok=True)

    for i, seed in enumerate(seeds.to_dict('records'), start=1):
        try:
            found, figure, diagnostics = evaluate_full_seed(
                seed, i, plot=not args.no_plots)
            events.extend(found)
            audit.append(dict(seed_utc=seed['start'],
                              seed_category=seed['category'],
                              detected_events=len(found),
                              figure=figure, error='', **diagnostics))
            print(f'[FULL {i}/{len(seeds)}] {seed["start"]}: '
                  f'{len(found)} candidate events', flush=True)
        except Exception as exc:
            print(f'[FULL ERROR] {seed["start"]}: {exc}', flush=True)
            audit.append(dict(seed_utc=seed['start'],
                              seed_category=seed['category'],
                              detected_events=0, figure='', error=str(exc)))
        save_progress()
        print(f'[SAVED] Analyzed {i}/{len(seeds)} windows; '
              f'current catalog: {output}', flush=True)
    # Avoid double-counting the same event in overlapping reader windows.
    unique = deduplicate(events)
    good = sum(not r['error'] for r in audit)
    print(f'COMPLETE: {len(seeds)} screening windows; '
          f'{len(events)} candidate event detections; '
          f'{len(unique)} after dedup; {good}/{len(seeds)} windows processed')
    print(f'Saved {output} and {audit_output}')
    if not good:
        raise SystemExit('No windows analyzed successfully: see audit CSV')


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--top',type=int,default=None,help='Windows to analyze; 0=all (full mode default: all; legacy default: 8)')
    p.add_argument('--mode',choices=('full','legacy'),default='full',
                   help='Full-cadence candidate analysis (default) or original 4cyc benchmark')
    p.add_argument('--reference',choices=('local','parker-spc'),default='local')
    p.add_argument('--self-test',action='store_true')
    p.add_argument('--no-plots',action='store_true')
    p.add_argument('--windows-csv',type=Path,default=None,
                   help='Reader candidate CSV in full mode; legacy mode defaults to built-in 20 windows')
    a=p.parse_args()
    if a.self_test:
        self_test()
        full_self_test()
        return
    if a.top is not None and a.top < 0:
        p.error('--top must be nonnegative')
    if a.top is None:
        a.top = 0 if a.mode == 'full' else 8
    if a.mode == 'full':
        run_full_analysis(a, p)
        return
    # Default candidates are stored here; a separate file is optional.
    if a.windows_csv is None:
        seeds=default_windows()
    else:
        if not a.windows_csv.is_file():
            p.error(f'Windows CSV does not exist: {a.windows_csv}')
        seeds=pd.read_csv(a.windows_csv)
        required={'start','end','category'}
        if not required.issubset(seeds.columns):
            p.error('Windows CSV must contain start,end,category columns')
    seeds['priority']=seeds.category.map({'strong':0,'ambiguous':1,'control':2}).fillna(3)
    seeds=seeds.sort_values(['priority','start'])
    if a.top:
        seeds=seeds.head(a.top)
    records=[];audit=[]
    for i,seed in enumerate(seeds.to_dict('records'),start=1):
        try:
            found,image=evaluate_seed(seed,i,reference=a.reference,plot=not a.no_plots)
            records.extend(found)
            audit.append(dict(seed_utc=seed['start'],seed_category=seed['category'],
                              detections_before_dedup=len(found),figure=image,error=''))
            print(f'[{i}/{len(seeds)}] {seed["start"]}: {len(found)} preliminary intervals',flush=True)
        except Exception as exc:
            print(f'[ERROR] {seed["start"]}: {exc}',flush=True)
            audit.append(dict(seed_utc=seed['start'],seed_category=seed['category'],
                              detections_before_dedup=0,figure='',error=str(exc)))
    unique=deduplicate(records)
    DATA.mkdir(exist_ok=True)
    name=a.reference.replace('-','_')
    pd.DataFrame(audit).to_csv(DATA/f'corrected_{name}_audit.csv',index=False)
    pd.DataFrame(unique,columns=[
        'seed_utc','seed_category','event_start_utc','event_end_utc','observed_span_seconds',
        'duration_seconds','left_censored','right_censored','maximum_deflection_deg',
        'median_abs_Bmag_deviation_pct','max_abs_Bmag_deviation_pct',
        'radial_reversal','reference','status','duplicate_source_count','source_seed_utcs'
    ]).to_csv(DATA/f'corrected_{name}_events.csv',index=False)
    censored=sum(bool(e['left_censored'] or e['right_censored']) for e in unique)
    print(f'COMPLETE: {len(seeds)} windows; {len(records)} raw detections; '
          f'{len(unique)} deduplicated; {censored} incomplete events')
    print(f'Files: data/corrected_{name}_events.csv and *_audit.csv')
    if not any(not row['error'] for row in audit):
        raise SystemExit('All windows failed. Inspect audit CSV; no result validated.')

if __name__=='__main__':
    main()
