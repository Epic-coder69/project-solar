# Independent holdout: first completed results (2026-10-25)

**Scope.** Frozen November 2018 Encounter 1 holdout selected with random seed
**41801** from the Martinović et al. visually annotated switchbacks and
catalog-clear unlabeled intervals. Selection excluded the original survey seed
exposures; the [manifest](holdout_manifest.csv) and
[settings](holdout_settings.json) preserve selection provenance.

## Data quality and coverage

- **20** distinct held-out intervals: **10 published positives** and
  **10 catalog-clear unlabeled** 30-minute periods.
- Both detection methods processed **all 20 intervals successfully**.
  All **40 per-method audit rows** have `status=ok`; the minimum reported
  fraction of valid one-second angular measurements is **0.9872**.
- The audit's high data coverage is encouraging; it **does not independently
  validate instrument quality flags or plasma-derived switchback physics**.

## Detector outcome on published examples

| Criterion | Local magnetic reference | Approximate Parker-SPC reference |
|---|---:|---:|
| Any candidate temporally overlaps the reference SB | 10/10 | 10/10 |
| IoU >= 0.25 with published SB boundaries | 2/10 | 2/10 |
| IoU >= 0.50 with published SB boundaries | 2/10 | 1/10 |

The **two** examples passing IoU >= 0.25 for both methods are
**H003 / E1-104** and **H006 / E1-129**:

| Example | Published duration | Local best IoU | Parker-SPC best IoU |
|---|---:|---:|---:|
| H003 / E1-104 | 868 s | 0.5474 | 0.2613 |
| H006 / E1-129 | 1428 s | 0.6635 | 0.5340 |

The common **10/10 overlapping-candidate** result must **not** be read as
validated event recovery: candidate intervals can be far longer than the
reference SB. Some are also window-censored. For instance, H001 has censored
intervals in both detector catalogs; H019's local interval is right-censored.
All three should be marked boundary-incomplete when interpreting duration.

Two positive examples (**H009/E1-156** and **H015/E1-394**) each last only
~**3.50 seconds**, shorter than the current detector's 10-second minimum
event duration. They are necessarily poor examples for judging strict
one-event-to-one-event segmentation using these fixed detector settings.

The H003 and H006 raw FIELDS figures show strong changes in the RTN
components and appreciable variations in scalar |B|; they should be
reviewed with available SWEAP plasma measurements before assigning physical
switchback status. H009 also shows complex field activity surrounding its
very short cataloged event. These are **qualitative figure observations, not
independent switchback classifications**.

## Outcomes on catalog-clear unlabeled intervals

Both methods flagged at least one magnetic rotation in the **same 7/10**
catalog-clear intervals: H005, H008, H010, H014, H017, H018, H020.
The remaining three intervals (H004, H007, H012) were not flagged by either.

These are **NOT validated false positives**. The Martinović catalog is
visually selected, not an exhaustive declaration that every unlisted
interval is event-free. For example, the H005 and H017 raw magnetic profiles
contain component variation and magnitude changes; automatic angular
threshold crossings there deserve physical review, not immediate rejection.

## Critical limitations

1. **Neither detector wins this test.** The 2/10 tie at IoU>=0.25 is
   particularly relevant to event-boundary analysis; at IoU>=0.50,
   the difference is only one example, with high uncertainty.
2. **Physical classification is outstanding.** The
   [blinded review sheet](holdout_blind_review.csv) still has **zero
   manual classifications**, so sensitivity/specificity/false-positive
   claims are unsupported.
3. **Catalog/event-definition mismatch matters.** The published catalog
   includes short intervals and joint magnetic/plasma visual selection.
   Our algorithms are magnetic-only with 10-second minimum and
   20-second gap merging.
4. **Bias and data availability:** The holdout avoids original seed windows
   but still samples positives from a curated published catalog and
   catalog-clear times from its coverage span. It is not a random draw
   from all PSP measurements.
5. **Censored detections:** The present best-IoU summary includes some
   edge-censored candidate spans, whose complete boundaries are unknown.

## Next decision and required review

Before choosing a detector, have an independent reviewer inspect the
[raw FIELDS review figures](../figures/holdout_review/) **without opening
the stratum manifest or algorithm results**. Capture reviewer notes in
`results/holdout_blind_review.csv`, preserving
`switchback_like` / `other_structure` / `indeterminate` distinctions.
Where possible add SWEAP plasma-speed and Alfvénic behavior to support
physical classification. Report boundary agreement only for complete,
quality-screened cases and separately characterize short events.

**Conclusion:** The two reference approaches both screen for strong magnetic
activity but neither achieves consistently accurate segmentation of the
independently published SB intervals. Their relative ranking remains
inconclusive. A detector has **not** been selected or frozen.

**Inputs:**
[`holdout_manifest.csv`](holdout_manifest.csv),
[`holdout_audit.csv`](holdout_audit.csv),
[`holdout_detections.csv`](holdout_detections.csv), and
[`holdout_summary.csv`](holdout_summary.csv).
