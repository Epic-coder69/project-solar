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

# One-to-one greedy matching mistakes split/merged intervals for exclusive detections.
# Compare groups of overlapping intervals from the same successfully evaluated seed windows.
parent = list(range(len(records['local']) + len(records['parker_spc'])))

def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

def connect(x, y):
    parent[find(x)] = find(y)

edges = []
for i, (ls, le) in enumerate(records['local']):
    for j, (ps, pe) in enumerate(records['parker_spc']):
        seconds = max(0, (min(le, pe) - max(ls, ps)).total_seconds() + 1)
        shorter = min((le - ls).total_seconds() + 1, (pe - ps).total_seconds() + 1)
        fraction = seconds / shorter if shorter else 0
        if fraction >= 0.5:
            connect(i, len(records['local']) + j)
            edges.append(dict(local_start=ls.isoformat(), parker_spc_start=ps.isoformat(),
                              overlap_seconds=round(seconds, 1),
                              shorter_overlap_fraction=round(fraction, 4)))

components = {}
for method in ('local', 'parker_spc'):
    offset = 0 if method == 'local' else len(records['local'])
    for index, interval in enumerate(records[method]):
        components.setdefault(find(offset + index), []).append((method, interval))

def pattern(n_local, n_spc):
    if not n_local: return 'parker_spc_only'
    if not n_spc: return 'local_only'
    if n_local == 1 and n_spc == 1: return 'one_to_one'
    if n_local > 1 and n_spc == 1: return 'local_splits'
    if n_local == 1 and n_spc > 1: return 'parker_spc_splits'
    return 'many_to_many'

def label(intervals):
    return '; '.join(s.strftime('%H:%M:%S') + '–' + e.strftime('%H:%M:%S')
                     for s, e in sorted(intervals))

groups = []
for items in sorted(components.values(), key=lambda es: min(iv[0] for _, iv in es)):
    local = sorted([iv for method, iv in items if method == 'local'])
    spc = sorted([iv for method, iv in items if method == 'parker_spc'])
    begin = min(iv[0] for _, iv in items)
    end = max(iv[1] for _, iv in items)
    groups.append(dict(group_id=len(groups) + 1, date=begin.date().isoformat(),
                       local_count=len(local), parker_spc_count=len(spc),
                       relationship=pattern(len(local), len(spc)),
                       group_start=begin.isoformat(), group_end=end.isoformat(),
                       local_intervals=label(local), parker_spc_intervals=label(spc)))

pd.DataFrame(groups, columns=['group_id', 'date', 'local_count', 'parker_spc_count',
                              'relationship', 'group_start', 'group_end',
                              'local_intervals', 'parker_spc_intervals']).to_csv(
                                  DATA/'reference_overlap_groups.csv', index=False)
pd.DataFrame(edges, columns=['local_start', 'parker_spc_start', 'overlap_seconds',
                             'shorter_overlap_fraction']).to_csv(
                                 DATA/'reference_overlap_pairs.csv', index=False)

counts = pd.Series([g['relationship'] for g in groups]).value_counts()
print('Shared windows:', len(shared))
print('Local/SPC detections:', len(records['local']), '/', len(records['parker_spc']))
print('Overlap-connected groups:', len(groups))
for name in ['one_to_one', 'local_splits', 'parker_spc_splits', 'many_to_many',
             'local_only', 'parker_spc_only']:
    print(name, int(counts.get(name, 0)))
print('IMPORTANT: shared time is not independent validation of physical switchbacks.')
print('Saved:', DATA/'reference_overlap_groups.csv', DATA/'reference_overlap_pairs.csv')
