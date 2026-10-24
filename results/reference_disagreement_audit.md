# Detector-reference disagreement audit — November 2018

**Data:** `results/corrected_local_events.csv`, `results/corrected_parker_spc_events.csv`, and their audit catalogs.

The comparison includes only the **three successfully evaluated observation windows shared by both methods** (November 1, 5 and 7, 2018). Both catalogs contain **10 candidate event intervals** in those shared windows.

## Temporal-overlap grouping

Instead of greedily requiring one-to-one event matches, detections are connected when they overlap in time by at least **50% of the shorter interval**. Connected sets form groups, which reveal different event segmentation.

| UTC date | Local intervals | Parker-SPC intervals | Relation | Approximate group period (UTC) |
|---|---:|---:|---|---|
| 2018-11-01 | 2 | 1 | Local splits | 18:20–19:17 |
| 2018-11-05 | 2 | 1 | Local splits | 13:34–13:41 |
| 2018-11-05 | 1 | 1 | One-to-one | 13:42 |
| 2018-11-05 | 1 | 2 | Parker-SPC splits | 13:44–13:58 |
| 2018-11-05 | 1 | 2 | Parker-SPC splits | 13:59–14:10 |
| 2018-11-07 | 1 | 1 | One-to-one | 03:18–03:57 |
| 2018-11-07 | 1 | 1 | One-to-one | 04:00–04:05 |
| 2018-11-07 | 1 | 1 | One-to-one | 04:07–04:09 |

**Summary:** eight connected groups: four one-to-one, two where the local detector splits a longer Parker-SPC interval, and two where Parker-SPC splits a longer local interval. At this overlap threshold, there are **zero completely unmatched groups**.

The older `corrected_reference_comparison.csv` reported eight one-to-one matches plus two `local_only` and two `parker_only` intervals. Those exclusive labels were artifacts of its greedy one-to-one matching constraint, because its supposedly exclusive intervals actually overlap detections in the other catalog.

## Physical interpretation and limitations

This is **agreement in the timing of magnetic activity, not independent validation of switchbacks**. Differences in event boundaries and whether deflections are split or merged are genuine algorithmic distinctions, and might change duration and rate statistics.

Neither algorithm is established as ground truth. The Parker-SPC mode is an exploratory approximation and is not the complete published Laker et al. SPAN-i-based background-reference method. Only three preselected windows are covered. No false-positive rate, missed-event rate, or population occurrence rate can be inferred.

**Next diagnostic test:** Overlay the magnetic deflection angle from both reference calculations on the same axes for each of the four segmentation-disagreement groups, mark the 30-degree boundary and 45-degree trigger, inspect local magnitude and data gaps, and compare with plasma velocities where available. Then repeat on independently chosen intervals.

**Reproducibility:** `compare_detectors.py` now produces `data/reference_overlap_groups.csv` and `data/reference_overlap_pairs.csv` from the local event catalogs, using a 50%-of-shorter overlap rule.
