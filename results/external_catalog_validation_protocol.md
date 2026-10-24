# Independent-catalog switchback validation protocol

## Purpose and provenance

This step tests temporal agreement with the **Martinović et al. (2021)** visually selected PSP Encounter 1–2 catalog (Zenodo DOI [10.5281/zenodo.4579821](https://doi.org/10.5281/zenodo.4579821)). The original file is `ApJ_SB_Database.h5` (published MD5 `5046a6964a2f17383c243892202ad5c5`). We do **not** include the third-party HDF5 in this Git repository.

The source catalog contains 1,074 visually selected events using magnetic-field rotations and bulk plasma velocity enhancements, not an exhaustive positive/negative truth set. Our selected November 2018 windows were chosen through magnetic screening, not independent random sampling. Accordingly, none of these summary fractions should be presented as *sensitivity, precision, specificity, false-positive rate, or unbiased occurrence rate*.

## Reproduce the analysis on the research Mac

Run from repository root after pulling latest changes:

```bash
git pull origin main
source .venv/bin/activate
python -m pip install h5py requests
python external_switchback_validation.py --self-test
python external_switchback_validation.py --inspect
```

The inspection step downloads the **original** Zenodo file into `data/external/` and prints its HDF5 keys, array shapes, dtypes, and first values. The downloadable HDF5 could not be inspected in the environment where this script was written. **Do not infer the event timestamp interpretation solely from numeric size.** Confirm the SB boundaries and timestamp units from the original file and its Zenodo ReadMe. Select the HDF5 table containing one row per published switchback.

Next run, replacing the illustrative key and field labels below with the **actual** labels from inspection:

```bash
python external_switchback_validation.py \
  --h5-key "ACTUAL_HDF5_DATASET_OR_GROUP" \
  --start-field "ACTUAL_SB_START_FIELD" \
  --end-field "ACTUAL_SB_END_FIELD" \
  --time-unit iso \
  --scope shared
```

Supported time interpretations are `iso`, `unix_s`, `unix_ms`, `unix_ns`, `mjd`, and `jd`. Do not use the default `iso` for numeric timestamps. When HDF5 storage is more complex, first extract a *verified* two-column event table to CSV, then use:

```bash
python external_switchback_validation.py \
  --catalog-csv data/external/verified_martinovic_events.csv \
  --start-field start_utc --end-field end_utc \
  --time-unit iso --scope shared
```

Inputs required locally: `data/survey_selected.csv` produced by `switchback_survey.py`, and detector outputs `results/corrected_{local,parker_spc}_{events,audit}.csv`.

## Evaluation decisions

- Default comparison includes only **shared, successfully analyzed seed windows** from the local and Parker-SPC reference algorithms.
- Reconstruct each seed's actual analysis interval using the `start` and `end` columns plus **40 minutes** on either side, matching `reference_detector_v2.py`. Merge overlapping windows.
- Only published reference events *fully contained* in the shared exposure qualify for the denominator. This excludes cases truncated by observation-window limits.
- A detection matches an external positive if there is a nonzero overlap equal to at least 50% of the **shorter** event interval. A maximum-cardinality **one-to-one** assignment avoids counting one broad detection as multiple recovered events.
- Save paired time intervals and overlap fractions for manual boundary review. A lenient shorter-overlap rule can match intervals with very different boundaries; report both reference/detector overlap fractions and test other thresholds.
- The default scope and time windows compare **candidate detection concordance with published positives only**. They do not measure undetected switchbacks outside the preselected windows or detector specificity.

## Outputs

The command writes in `results/`:

- `external_reference_eligible.csv`: published positives wholly inside the shared exposure
- `external_catalog_matches.csv`: matched positive–candidate pairs, timestamps and overlap fractions
- `external_catalog_summary.csv` and `.json`: counts and descriptive agreement fractions with the exact windows and matching rule

Before making scientific claims, verify that the parsed source has plausible UTC times and **the expected published row count**; compare matching figures against original spacecraft observations, and repeat on independent non-event intervals. The new matching code is infrastructure, **not yet an independent result** until the real Zenodo HDF5 dataset is parsed and analyzed.
