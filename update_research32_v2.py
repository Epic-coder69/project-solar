#!/usr/bin/env python3
"""Replace only research Section 3.2 with measured corrected benchmark results."""
import re
from pathlib import Path
import pandas as pd

paper=Path('The_Research.md')
local=Path('data/corrected_local_events.csv')
audit_path=Path('data/corrected_local_audit.csv')
if not paper.exists() or not local.exists() or not audit_path.exists():
    raise SystemExit('Missing The_Research.md or corrected local catalogs; run v2 first.')
text=paper.read_text(encoding='utf-8')
start=re.search(r'(?m)^###\s+3\.2\b[^\n]*$',text)
end=re.search(r'(?m)^###\s+3\.3\b[^\n]*$',text)
if start is None or end is None or end.start()<=start.start():
    raise SystemExit('Section 3.2 / 3.3 headings not found. Nothing changed.')
events=pd.read_csv(local)
audit=pd.read_csv(audit_path)
if audit.empty or not audit['error'].astype('string').fillna('').eq('').any():
    raise SystemExit('No successfully evaluated intervals. Paper unchanged.')
success=audit['error'].astype('string').fillna('').eq('').sum()
censored=(events.left_censored.fillna(False).astype(bool) | events.right_censored.fillna(False).astype(bool)) if len(events) else pd.Series(dtype=bool)
complete=events.loc[~censored] if len(events) else events
repeat=int((events.duplicate_source_count.fillna(1)-1).sum()) if len(events) else 0
median_dur=f'{complete.duration_seconds.median():.0f} seconds' if len(complete) and complete.duration_seconds.notna().any() else 'not established'
median_dev=f'{events.median_abs_Bmag_deviation_pct.median():.2f}%' if len(events) and events.median_abs_Bmag_deviation_pct.notna().any() else 'not established'
parker_file=Path('data/corrected_parker_spc_events.csv')
parker_line=('An exploratory SWEAP-SPC Parker-spiral-reference comparison was also run. '
             f'Its catalog contains **{len(pd.read_csv(parker_file))} unvalidated intervals**. '
             'It is NOT the full Laker et al. method.' if parker_file.exists() else
             'A SWEAP-SPC Parker-spiral comparison is available in code, but has not yet produced a verified catalog.')
figure_lines=[]
for d in ('20181107_0332','20181108_0944'):
    files=list(Path('figures').glob(f'corrected_local_strong_{d}_*.png'))
    if files:
        figure_lines.append(f'![Corrected local-background detection at {d}](figures/{files[0].name})')
figures='\n\n'.join(figure_lines) or '*Diagnostic images will be added after the corrected run.*'
section=f'''### 3.2 Switchback Identification and Characterization

#### 3.2.1 Detector and reference definitions

Our benchmark retains the angular detection thresholds used by Laker et al. (2022):
45° trigger, 30° boundary, 10-second minimum duration, and merging separated
segments when the intervening valid-data gap is less than 20 seconds.

For this **local-background benchmark**, the magnetic deflection angle is
calculated against a 181-minute rolling median of the magnetic RTN components:

$$
\\theta(t)=\\cos^{{-1}}\\left(\\frac{{\\mathbf{{B}}(t)\\cdot\\mathbf{{B}}_{{bg}}(t)}}{{|\\mathbf{{B}}(t)|\\,|\\mathbf{{B}}_{{bg}}(t)|}}\\right).
$$

The independent magnitude reference is the 181-minute rolling median of **scalar**
$|\\mathbf{{B}}|$, not the magnitude of the component-wise median vector.
Magnitude deviations are reported as percentages of this scalar reference.

This benchmark is **not an exact reproduction** of Laker et al.'s Parker-spiral
method. The optional `parker-spc` mode uses good-quality SWEAP SPC Level 3
primary-proton radial velocity, spacecraft HCI position, an assumed 2.5 solar-radius
source surface and local magnetic-sector polarity. It is exploratory: the published
method primarily uses SPAN-i measurements and a filtered velocity reference.

#### 3.2.2 Event catalog and data-quality safeguards

The detector uses one-second median magnetic samples; missing samples are not
interpolated. Each event records a peak angle, estimated start and end, radial
reversal, scalar-field-magnitude deviation, and reference method.

An event that begins or ends at the analysis-window boundary (or adjacent
to a missing angular sample) is **censored**: its complete duration is unknown.
Its `duration_seconds` is left blank, while `observed_span_seconds` is retained.
Duplicated detections from overlapping seed windows are consolidated in the
corrected catalog. All detections remain **unvalidated magnetic rotations**.

#### 3.2.3 Corrected preliminary results

- Selected windows successfully processed: **{success}**.
- Original seed-associated detections after global deduplication: **{len(events)}**.
- Duplicate observations consolidated: **{repeat}**.
- Boundary/gap-censored detections: **{int(censored.sum())}**.
- Median duration among only uncensored detected intervals: **{median_dur}**.
- Median absolute scalar $|B|$ deviation across detected intervals: **{median_dev}**.

These values describe a deliberately selected set of observation windows. They are
**not** population-wide occurrence rates, event-confirmation statistics, or an
estimate of false-positive rate.

{parker_line}

#### 3.2.4 Selected diagnostic figures

{figures}

#### 3.2.5 Remaining validation

A full replication must reproduce the published plasma-based Parker-spiral
reference (including SPAN-i primary measurements and smoothing), apply consistent
quality flags, test threshold sensitivity, independently review candidates and
non-candidates, and quantify false positives and missed events before scaling to
multi-year statistical results.

**Reference:** [Laker et al. (2022)](https://academic.oup.com/mnras/article/517/1/1001/6692877).

'''
backup=paper.with_suffix('.md.pre_v2.bak')
backup.write_text(text,encoding='utf-8')
paper.write_text(text[:start.start()]+section+text[end.start():],encoding='utf-8')
print('Updated Section 3.2 (backup:',backup,')')
print(f'Unvalidated unique intervals={len(events)}, censored={int(censored.sum())}; no population-level conclusions.')
