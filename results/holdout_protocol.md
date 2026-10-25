# Independent switchback holdout experiment

## Scope and scientific caution

The runner at [independent_holdout.py](../independent_holdout.py) selects intervals independently of either detector output. Its reference is the visually selected [Martinović et al. PSP E1–E2 catalog](https://doi.org/10.5281/zenodo.4579821), which contains 1,074 switchbacks in five-region records.

We sample two distinct groups in November 2018: (1) ten published positives selected randomly outside the original 20 observation windows, and (2) ten randomly timed **catalog-clear, UNLABELED** intervals within the Encounter-1 annotated date span. If the exclusion rules prevent reaching either count, the manifest records the smaller sample. Deterministic RNG seed: **41801**.

**No leakage:** The entire original survey selection plus 40-minute detector context on each side is excluded. The new evaluation spans cannot overlap each other. Published-positive sampling depends only on external published timestamps and original survey seed positions, never on our previous detections.

**Unlabeled control:** A 30-minute center interval with a 10-minute buffer around all known published event regions qualifies as catalog-clear. It may still contain an event the catalog did not annotate. Days need at least 1,000 valid 1-minute observations in the original coverage summary. Data coverage is finally checked at 1-second resolution when the detector runs.

## Run from the project directory on your Mac

    git pull --ff-only origin main
    source .venv/bin/activate
    python -m pip install h5py numpy pandas matplotlib requests cdflib
    python independent_holdout.py self-test
    python independent_holdout.py prepare
    python independent_holdout.py run --plots

The prepare command freezes a SHA-256 verified results/holdout_manifest.csv and results/holdout_settings.json. It refuses to resample an existing holdout to avoid post-result cherry picking. Both detector references use the unchanged v2 detection thresholds and the same target intervals. Each must have at least 80% valid 1-second angle samples in the evaluation region; unavailable data are marked as failures, not non-events.

## Saved evidence

- results/holdout_manifest.csv — holdout IDs, source strata and original published reference IDs. Keep separate from blind review.
- results/holdout_settings.json — sampling seed, provenance, and manifest SHA-256.
- results/holdout_audit.csv — each method's success, errors, and observed valid-data fraction.
- results/holdout_detections.csv — candidate bounds and peak angles.
- results/holdout_summary.csv — candidate counts and intersection-over-union (IoU) matches to published positives at thresholds 0.25 and 0.50, strictly on paired valid windows.
- results/holdout_blind_review.csv — **blank reviewer classifications** and figures without detector labels or knowledge of positive versus catalog-clear status.
- figures/holdout_review/H*.png — RTN magnetic-field components and magnitude, fixed-width central review region.

After completing the Mac run:

    git add results/holdout_*.csv results/holdout_*.json figures/holdout_review/
    git commit -m "Run independent held-out switchback detector comparison"
    git push origin main

## Reviewer and interpretation protocol

Before opening the manifest or detector outputs, inspect the raw figures. Enter one of **switchback_like**, **other_structure**, or **indeterminate** in the manual_classification column and record notes. The raw-FIELDS review is only a magnetic-structure classification: true physical switchbacks may require SWEAP plasma-speed/Alfvénic checks.

Published positives are independently selected examples, not an exhaustive truth set. Catalog-clear intervals are not confirmed negatives. The two groups have deliberate selection bias and cannot be combined to estimate mission-wide precision, sensitivity, specificity, false-positive rate, or occurrence rate. Compare the event-boundary IoU, availability of the SPC plasma reference, and descriptive detector flagging separately. Do not tune thresholds after inspecting these holdout results; a larger, multi-encounter expert-labeled dataset would still be needed for definitive validation.
