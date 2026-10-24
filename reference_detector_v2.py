#!/usr/bin/env python3
"""Project Solar: corrected FIELDS switchback-candidate benchmark.

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

Requires psp_fields.py and switchback_survey.py from the existing project.
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

from psp_fields import get_day
from switchback_survey import read_high

DATA = Path('data')
FIG = Path('figures')
SPC_CACHE = DATA / 'sweap_spc'
SPC_API = ('https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/sp_phys/'
           'datasets/PSP_SWP_SPC_L3I/orig_data')
TRIGGER_DEG, EDGE_DEG = 45., 30.
MIN_SECONDS, MERGE_SECONDS = 10, 20
SOLAR_RADIUS_KM = 695700.
OMEGA_SUN_RAD_S = 2 * np.pi / (25.38 * 86400.)

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


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--top',type=int,default=8,help='Seed windows to analyze; 0=all')
    p.add_argument('--reference',choices=('local','parker-spc'),default='local')
    p.add_argument('--self-test',action='store_true')
    p.add_argument('--no-plots',action='store_true')
    a=p.parse_args()
    if a.self_test:
        self_test();return
    if a.top<0:
        p.error('--top must be nonnegative')
    catalog=DATA/'survey_selected.csv'
    if not catalog.exists():
        p.error('Missing data/survey_selected.csv. Run existing survey first.')
    seeds=pd.read_csv(catalog)
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
