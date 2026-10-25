# Independent catalog cross-match: November 2018 — methodological audit

**Source:** [Martinović et al., visually selected PSP E1–E2 switchbacks](https://doi.org/10.5281/zenodo.4579821).
**Data:** `results/external_catalog_summary.json`, `external_reference_eligible.csv`,
`external_catalog_matches.csv`, and both `corrected_*` detection catalogs.
**Analysis recorded:** 2026-10-24. **Scope:** 3 shared, successfully evaluated seed windows,
not an unbiased spacecraft-wide sample.

## Observed results of the original matching criterion

The publisher's 1,074-switchback catalog contributes **17 central switchback
intervals** fully contained in the three analyzed windows
(2018-11-01 18:02:30–19:23:30; 2018-11-05 13:10:30–14:35:30;
2018-11-07 02:52:30–04:34:30 UTC).

| Metric | Local reference | Parker-SPC reference |
|---|---:|---:|
| Candidate intervals in scope | 10 | 10 |
| Published events matched using >=50% of *shorter* duration | 6 | 7 |
| Matched fraction of 17 published events | 35.3% | 41.2% |

These are **catalog overlap fractions**, not sensitivity, recall, precision,
specificity, or unbiased event rates. The windows were preselected by our
original magnetic screening, and the external catalog is a visually selected,
positive-only dataset.

## Serious issue with the original matching rule

Counting half the **shorter** interval as sufficient ignores excessive
duration of the longer event:

- Published event **E1-233** lasts ~3.50 seconds (2018-11-05 14:01:25–14:01:28 UTC).
  Both detectors match it to a broad ~9–11-minute candidate even though the
  common interval is only about **0.5–0.6%** of the detector interval.
- Published event **E1-338** lasts ~58.31 seconds (2018-11-07 03:20:14–03:21:12 UTC).
  Both match it to a ~38-minute candidate, covering only **~2.5–2.6%**
  of the detector interval.
- The published E1-340 switchback (2018-11-07 03:33:27–03:55:44 UTC)
  substantially overlaps the much longer detector candidate in that interval.
  More balanced matching can prefer this event instead of E1-338.

Thus the apparent 7-to-6 advantage may reward **overly broad event boundaries**
rather than accurate event identification.

## Threshold sensitivity on the same 17 published intervals

Recomputed from the saved positive and detection catalogs using maximum-cardinality
one-to-one temporal matches. Intersection-over-union (IoU) is
`duration(intersection) / duration(union)`.

| Event matching rule | Local matches | Parker-SPC matches |
|---|---:|---:|
| >=50% of shorter interval (original) | 6/17 | 7/17 |
| Reference >=50% AND detector >=25% | 5/17 | 5/17 |
| IoU >=0.10 | 5/17 | 7/17 |
| **IoU >=0.25** | **5/17** | **5/17** |
| IoU >=0.50 | 4/17 | 3/17 |
| Overlap >=30 sec; reference >=50%, detector >=25% | 4/17 | 4/17 |

**There is no robust winner.** Small differences reverse or disappear with
reasonable variations in how event boundaries are evaluated.
The new `external_match_sensitivity.py` reproduces this analysis and exports
the matched pairs for each rule.

## Duration and definition limitations

The 17 reference intervals have a wide range of durations.
**Seven are shorter than 30 seconds** and **three are shorter than 10 seconds**.
Our current detector imposes `MIN_SECONDS = 10` and merges qualifying
threshold crossings with `MERGE_SECONDS = 20` (see
`reference_detector_v2.py`), so it cannot isolate every short reference
interval using its present definitions. The published catalog also relies on
velocity behavior, whereas the current detector is magnetic-only.

These are different **event definitions** and segmentation choices, not merely
technical mistakes. Do not tune the thresholds repeatedly to these same
17 events and then report performance on them as an independent test.

## Next scientific decision

Neither method should be declared superior yet.

1. Inspect the unmatched and broad-match intervals visually, comparing the
   published central SB boundaries, FIELDS magnitude and field direction,
   and plasma velocity where available.
2. Build an expanded, independently selected holdout set containing published
   positives and human-reviewed non-events, spanning multiple encounters,
   solar distances and data-quality conditions.
3. Report results across a prespecified range of IoU and duration thresholds;
   separate *identification of magnetic activity* from *event-boundary accuracy*
   and characterize the prevalence of short events.
4. Only then select/freeze a detector for the population survey.

**Research status:** an independent positive-catalog comparison has been run,
but has not established detector accuracy or a final detection method.


## Expanded shared-window comparison — October 25, 2026

The subsequent run analyzed **16 of the 20 selected survey seed windows**
with both methods. Four Parker-SPC seed windows had insufficient reference
coverage (2018-11-25 21:29, 2018-11-26 20:25, 2018-11-26 22:20, and
2018-11-27 04:33 UTC) and were **excluded from the common-window comparison**.
The 16 shared seeds merge to **14 exposure intervals**, containing **27**
published central switchback intervals fully within their boundaries.

The local method produced **48** candidate intervals in the shared exposures;
the Parker-SPC approximation produced **44**. Results are reproducible from
[`external_catalog_summary.json`](external_catalog_summary.json),
[`external_match_sensitivity.csv`](external_match_sensitivity.csv), and the
individual [matched pairs](external_match_sensitivity_pairs.csv).

| Matching rule | Local reference | Parker-SPC approximation |
|---|---:|---:|
| At least 50% overlap of shorter event | 11/27 | 12/27 |
| At least 50% reference and 25% detector overlap | 8/27 | 8/27 |
| IoU >= 0.10 | 10/27 | 12/27 |
| **IoU >= 0.25** | **9/27** | **9/27** |
| IoU >= 0.50 | 7/27 | 5/27 |
| >= 30 seconds AND balanced overlap | 7/27 | 7/27 |

The difference remains **highly dependent on the definition of an event
match**, with no defensible winner. These counts measure overlap with
**visually selected positives in preselected windows**; they cannot estimate
mission-wide sensitivity, precision, specificity, or false-positive rate.

### Illustrative boundary agreement

Published event **E1-457** on November 8 has strong alignment with both
detectors: local IoU about **0.944**, Parker-SPC IoU about **0.947**.
Other detections have implausibly broad durations. In particular, the
published few-second **E1-233** and one-minute **E1-338** events can be counted
by the shorter-event overlap rule even when less than **3%** of the detector's
long interval overlaps the published event. This is why a balanced metric is
essential.

### Implications for the next experiment

1. Keep both methods as unvalidated **candidate screeners**, not competing
   established classifiers.
2. Compare published intervals with no detected events, and independently
   review representative negative/control periods and problematic boundaries.
3. Investigate four unusable Parker-SPC windows: the velocity reference was
   missing, not proof that switchbacks were absent.
4. Plan evaluation on **independently sampled observation intervals**, not
   only the windows prioritized by the initial detector, before freezing
   an algorithm for population-rate studies.
