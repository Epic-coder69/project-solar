#!/usr/bin/env python3
"""Compare reference algorithms on successfully evaluated *same* survey seeds.
No ground-truth classifications are inferred from agreement.
"""
from pathlib import Path
import pandas as pd

DATA = Path('data')
paths={'local':DATA/'corrected_local_events.csv',
       'parker_spc':DATA/'corrected_parker_spc_events.csv'}
audits={'local':DATA/'corrected_local_audit.csv',
        'parker_spc':DATA/'corrected_parker_spc_audit.csv'}
for p in list(paths.values())+list(audits.values()):
    if not p.exists():
        raise SystemExit(f'Missing {p}. Run the corresponding reference first.')

valid_seeds=[]
for method,p in audits.items():
    df=pd.read_csv(p)
    good=df.loc[df.error.astype('string').fillna('').eq(''),'seed_utc']
    valid_seeds.append(set(good))
shared=set.intersection(*valid_seeds)
if not shared:
    raise SystemExit('No jointly successful seed windows, so no fair comparison is possible.')
print('Shared successfully evaluated seed windows:',len(shared))

records={}
for method,p in paths.items():
    df=pd.read_csv(p)
    if df.empty:
        records[method]=[]
        continue
    def selected(row):
        sources=str(row.get('source_seed_utcs',row.seed_utc)).split(';')
        return bool(shared.intersection(sources))
    df=df[df.apply(selected,axis=1)]
    records[method]=[(pd.Timestamp(r.event_start_utc),pd.Timestamp(r.event_end_utc))
                     for r in df.itertuples(index=False)]

pairs=[];used=set()
for start,end in records['local']:
    best_score=0.;best_j=None
    for j,(other_start,other_end) in enumerate(records['parker_spc']):
        if j in used:continue
        overlap=(min(end,other_end)-max(start,other_start)).total_seconds()+1
        shorter=min((end-start).total_seconds()+1,(other_end-other_start).total_seconds()+1)
        score=max(0.,overlap)/shorter if shorter>0 else 0.
        if score>=.5 and score>best_score:
            best_j=j;best_score=score
    if best_j is None:
        pairs.append(dict(local_start=start.isoformat(),parker_start='',result='local_only'))
    else:
        used.add(best_j)
        pairs.append(dict(local_start=start.isoformat(),
                          parker_start=records['parker_spc'][best_j][0].isoformat(),
                          result='both'))
for j,(start,end) in enumerate(records['parker_spc']):
    if j not in used:
        pairs.append(dict(local_start='',parker_start=start.isoformat(),result='parker_only'))
out=DATA/'corrected_reference_comparison.csv'
pd.DataFrame(pairs,columns=['local_start','parker_start','result']).to_csv(out,index=False)
counts=pd.Series([x['result'] for x in pairs]).value_counts()
print(f'Compared local events: {len(records["local"])}; SPC Parker events: {len(records["parker_spc"])}')
for key in ('both','local_only','parker_only'):
    print(f'{key}: {int(counts.get(key,0))}')
print('Agreement only—NOT sensitivity, specificity, or verified switchbacks.')
print('Saved:',out)
