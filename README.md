# Project Solar

**Parker Solar Probe | Magnetic switchbacks, energetic particles, and plasma**

Project Solar investigates how magnetic-field structures identified in
Parker Solar Probe **FIELDS** data relate to energetic particles measured by
**IS⊙IS** and plasma conditions measured by **SWEAP**. The goal is a
reproducible, scientifically careful cross-instrument study — not merely an
automated switchback detector.

> **Start here:** run commands from the repository root using
> `python project.py`. This avoids path changes and the special
> character in `code/IS⊙IS/`.

## Current research status

Our November–December 2018 FIELDS catalog contains **599 automatically
selected magnetic switchback candidates**, not 599 independently confirmed
physical events. Our last completed IS⊙IS timing-coverage audit found:

| Observation tier | Candidates | What it actually supports |
| --- | ---: | --- |
| Fine EPI-Hi 1–10-second timing with event and ±10-second flanks | **295/599** | Candidate-scale timing study, pending particle quality/count checks |
| Additional HET proton/electron hourly-background coverage | **302/599** | Broader solar energetic-particle conditions; **not** switchback-resolved |
| No matching measurements in those audited HET channels | **2/599** | Missing, **not** zero particle activity |

For the original LET1 proton channel, **one** of the additional candidates
has partial coarse context instead of hourly-only context. The original
combined LET1 heavy-ion Z=6–28 channel has **597 hourly** observations and
**no established identical event-resolved replacement**.

**Important:** 597/599 (99.7%) represents nominal **at-least-hourly time
coverage**, *not* 99.7% high-resolution coverage or an independently
verified energetic-particle detection rate. Instrument quality flags,
passband compatibility, sparse-count statistics, viewing geometry,
selection bias and physical associations still need validation.

## Quick start

Requires Python 3.10+ and the dependencies in the tracked requirements file:

```bash
git clone https://github.com/Epic-coder69/project-solar.git
cd project-solar
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r code/FIELDS/external_libraries_requirements.txt
python project.py status
python project.py test
```

If you already have the repository and local NASA data, **do not reclone**.
Update only tracked code and documentation:

```bash
cd ~/Desktop/project-solar
git pull --ff-only origin main
source .venv/bin/activate
python project.py status
```

On macOS, to open the existing folder in VS Code without installing its
shell command: `open -a "Visual Studio Code" .`

### Unified IS⊙IS interfaces

There are now **exactly two IS⊙IS Python source files**, not wrappers:

- `code/IS⊙IS/isois_reader.py` contains both original L2 summary CDF
  decoding and EPI-Hi native-cadence CDF decoding, including their tests.
- `code/IS⊙IS/isois_analyzer.py` contains the original hourly survey
  analysis, fine-event coverage audit and multiresolution coverage audit,
  including their tests.

The five former `particle_*.py` scripts and the temporary interface-test
wrapper were merged and then removed. This reduces the IS⊙IS source folder to
**one actual reader and one actual analyzer**; the functions are in those two
files, not delegated to hidden subprocess backends.

- `isois-reader` **ingests** official CDF files, preserves the original
  detector, species, energy bins, native integration times, count/flux values,
  quality flags and source metadata. The `epihi` and `summary` modes are
  intentionally separate products behind the **same reader interface**.
- `isois-analyzer` **analyzes existing local measurements**: event-resolved
  coverage (`fine-coverage`), same-channel multi-resolution coverage
  (`coverage`) and the exploratory hourly survey (`survey`). The
  results keep the original CSV/JSON paths and schemas.

```bash
# Run all eight underlying offline instrument tests (no NASA downloads)
python project.py test

# Inspect the public interfaces
python project.py isois-reader --help
python project.py isois-reader list-products
python project.py isois-analyzer --help

# Read EPI-Hi observations at their NATIVE resolutions
python project.py isois-reader epihi \
  --events-csv data/full_analysis_nov_dec_2018_candidates_strong_events.csv \
  --products let1-rates3600 het-rates3600 \
  --download --resume

# Analyze the SAVED measurements, without redownloading any files
python project.py isois-analyzer coverage
python project.py isois-analyzer fine-coverage

# View the SWEAP products (separate instrument suite)
python project.py sweap --help
python project.py sweap --list-products
```

Legacy aliases **continue to work unchanged**, including
`python project.py highres ...`, `python project.py particles ...`,
`python project.py coverage`, `python project.py fine-coverage ...` and
`python project.py particle-analysis ...`. The former standalone `particle_*.py` script filenames are retired;
use `project.py` or the two new files instead. No cached NASA CDFs, NPZs,
summary CSVs or existing analysis results are moved, rewritten or deleted.

**Scope:** this is a structural refactor, **not** a new particle-physics
pipeline. The existing coverage reports measure nominal timestamp overlap;
they do *not* newly interpret instrument quality flags or establish Poisson
significance, energy-passband equivalence or particle acceleration. The
previous 295/599 fine and 597/599 at-least-hourly figures remain the last
verified user-run results until the analysis is rerun.

### Five-minute EPI-Hi heavy-ion rate recovery

The archived `let1-rates300` and `het-rates300` products contain
**sectored heavy-ion count rates** (CNO, Ne–Si, Fe-group), not `*_Flux`
arrays. Earlier reader code rejected these files because it incorrectly
required a differential-flux array. The corrected reader retains 300-second
integration boundaries, energy/sector axes, count-rate values, numeric
validity masks and source metadata, classifying these products as
`sectored_heavy_ion_count_rates_only`.

The analyzer separately reports `LET1_R1A_*`, `LET1_R26A_*` and
`HET_A_*` heavy-ion groups, **not** as the old Z=6–28 aggregate or as
interchangeable calibrated flux. Nominal timestamp overlap is not a
quality-flag or Poisson-significance check.

Reprocess already-cached NASA source CDFs **without downloading again**:

```bash
git pull --ff-only origin main
python project.py isois-reader --self-test
python project.py isois-analyzer --self-test
python project.py isois-reader epihi \
  --events-csv data/full_analysis_nov_dec_2018_candidates_strong_events.csv \
  --products let1-rates300 het-rates300 --offline --resume
python project.py isois-analyzer coverage
```

The 295 candidates with 1–10-second timing remain the same primary
high-resolution sample. Five-minute heavy-ion count rates are separate
composition/time-context measurements; they do not reconstruct missing
1–10-second proton or electron observations.

### Graph IS⊙IS measurements during FIELDS switchbacks

The unified analyzer now creates **native-value PNG plots** (no additional
Python scripts) with event start/end shaded. Each panel shows an independent
particle channel and a **single source product**. Horizontal line segments
cover the *actual integration interval* for each rate or differential-flux
value. Invalid/fill bins and missing data remain gaps; no interpolation,
cross-product averaging, flag interpretation or acceleration inference occurs.

From your activated project environment (plots require optional Matplotlib):

```bash
python -m pip install matplotlib
python project.py isois-analyzer --self-test
python project.py isois-analyzer plot-events --max-events 4
```

The first four plots prioritize switchbacks with fine HET proton coverage, if
the saved coverage report is present. PNGs and a provenance JSON summary are
written to `data/isois/plots/`. No NASA files are fetched, and the source
data in `data/` are not overwritten.

To plot a particular FIELDS switchback or different source channels:

```bash
python project.py isois-analyzer plot-events --list-events
python project.py isois-analyzer plot-events --event-id EVENT_ID --channels LET1_A_H HET_A_H HET_A_Electrons LET1_R26A_CNO
python project.py isois-analyzer plot-events --max-events 3 --energy-index 1 --sector-index 0
```

Energy and sector indices select one source-native bin; the analyzer never
silently sums unrelated energy bins or telescope directions. Hourly data are
labelled as such and must not be interpreted as switchback-scale variation.
See `python project.py isois-analyzer plot-events --help`.

## SWEAP reader (one file)

All SWEAP discovery, acquisition, single-product inspection and preliminary
eight-family survey functions now live in **one source file**:
[`code/SWEAP/sweap_data_reader.py`](code/SWEAP/sweap_data_reader.py).
There is no separate `sweap_family_survey.py` module. `--survey` is a
**reader option for preliminary data inspection**, not the scientific
plasma-analysis stage. A future `sweap_data_analysis.py` would deal with
quality-screened event measurements, cross-instrument comparisons, and
the 599 magnetic candidates.

The reader supports an **eight-family / thirteen-product** NASA SWEAP
inventory (SPC cup spectra and moments, SPAN-I proton/alpha distributions and
partial moments, SPAN-Ae/B full 3D electrons and energy spectra, and electron
pitch-angle distributions). Start by checking file availability with **no
downloads**:

```bash
python project.py sweap --self-test
python project.py sweap --survey --date 2018-11-06 --catalog-only
```

Then run a cautious first look at the actual multidimensional CDF data:

```bash
python project.py sweap --survey --date 2018-11-06 \
  --minutes 8 --max-file-mib 50 --max-total-mib 250 \
  --json data/sweap/survey_20181106.json
```

The first November 6 survey found **8 inspected products**, **5 original
daily files skipped by the 50 MiB cap**, and **0 processing errors**. Those
eight products represent **4 of the 8 measurement families**; skipped is
not missing. To investigate just the four oversized families, the reader
can request **NASA-generated eight-minute CDF subsets** instead of
large full-day CDFs:

```bash
python project.py sweap --survey --date 2018-11-06 \
  --families cup-charge-flux,proton-3d,alpha-3d,electron-3d \
  --subset-on-limit --minutes 8 \
  --max-file-mib 50 --max-total-mib 250 \
  --json data/sweap/survey_large_families_20181106.json
```

This fallback is opt-in; original daily files, including any already
downloaded files, are not overwritten. Subsets live under
`data/sweap/_subsets/`. If NASA cannot provide a subset, or the
subset exceeds the same download budget, the report records that condition
instead of silently omitting the product.

**Measured first-pass SWEAP survey outcome (2018-11-06):** The original
bounded run inspected **8/13 products** and skipped 5 oversized daily CDFs,
with 0 errors. The follow-up bounded NASA time-subset run inspected **all 5
remaining products**, with 0 errors or size-limit skips. Thus **all 13
product variants across the eight defined families have now been inspected
against real NASA CDFs**, although this does not establish scientifically
valid measurements for all 599 events.

| Product family | Variants | Records per tested 8-minute window |
| --- | --- | ---: |
| SPC L2 charge-flux spectra | 1 | 550 |
| SPC L3 proton moments | 1 | 549 |
| SPAN-I L2 proton flux (3D) | 1 | 17 |
| SPAN-I L2 alpha flux (3D, proton-contaminated) | 1 | 8 |
| SPAN-I L3 partial moments | Proton / alpha | 17 / 8 |
| SPAN-E L2 3D distributions | SPA / SPB | 17 / 17 |
| SPAN-E L2 summed energy spectra | SPA / SPB | 275 / 275 |
| SPAN-E L3 pitch-angle distributions | Combined / SPA / SPB | 17 / 17 / 17 |

**Interpretation:** These counts show EPOCH record availability in the chosen
test windows, not actual integration duration, uninterrupted valid coverage,
quality-flag pass rates, or co-observation of a FIELDS switchback. The
machine-local full survey reports are
`data/sweap/survey_20181106.json` and
`data/sweap/survey_large_families_20181106.json` (ignored by Git). Next:
review variable axes, quality flags, valid values and timestamp spacing,
then audit the specific candidate intervals without extrapolating one
eight-minute window to the entire mission.


This writes product-by-product variable names, units, array dimensions,
sampled numeric values, measured EPOCH spacing, and quality-flag variable
names. Missing days, unavailable products, oversized downloads, and errors
are reported explicitly rather than interpreted as zeros. It is a **survey**,
not a validated quality/coverage audit of the 599 switchback candidates.

## SWEAP: maximal native data capture at FIELDS switchbacks

The single `code/SWEAP/sweap_data_reader.py` now supports
`--extract-events` and `--audit-events`. It preserves every time-alignable
native variable, full instrument multidimensional distributions, original
fill values, uncertainty and quality flags, science and support axes, full
available CDF variable attributes and global metadata. It operates on the
exact half-open FIELDS `event_start_utc` –
`event_end_exclusive_utc` interval, never interpolating or averaging.

By default, the saved window spans the magnetic event **plus** 10 seconds
before and after and an additional **120 seconds of coarse-support context
on each side**. Those phases are marked separately in the metadata.
This extra context helps preserve low-cadence ~28–60-second spectra.
Neither the outer support nor an EPOCH inside the event is evidence of
actual measurement-integration overlap. The SPC L2 reader also derives an
additional `derived_SPC_step_UTC_ns` array from the original `Epoch` and
`MEASUREMENT_TIME` (nanosecond offsets per voltage step), without changing
either original variable. These are per-step timestamps, not exposure lengths.

Small source CDFs are cached locally. Large daily files are requested as
NASA-generated native-resolution time subsets; neighboring event requests
are grouped and reused when possible, with an individual-event retry
if an oversized grouped subset is unavailable. Download budgets are
configurable. There is **no default per-variable array size cap**
(`--max-variable-mib 0`), but sources may still be limited by the
per-file and total-download budgets, missing NASA coverage, or CDF
time-axis ambiguities. Every omission is documented.

**First real-data validation: just one event**

```bash
cd ~/Desktop/project-solar
git pull --ff-only origin main
python project.py sweap --self-test
python project.py sweap --extract-events \
  --events-csv data/full_analysis_nov_dec_2018_candidates_strong_events.csv \
  --max-events 1 \
  --before-seconds 10 --after-seconds 10 \
  --support-seconds 120 \
  --max-file-mib 250 --max-total-mib 1200 \
  --event-out data/sweap/events_v2
python project.py sweap --audit-events --bundle-audit --event-out data/sweap/events_v2
```

The output is isolated from the earlier prototype archives under
`data/sweap/events_v2/` and remains Git-ignored. For each product/event/UTC
day there is a losslessly sampled `.npz` archive (compressed, not resampled)
with an adjacent `.json` containing the source reference, native dimensions,
full CDF variable attributes, missing/invalid counts, and exact EPOCH phase
counts. The output also includes `event_manifest.csv`, `summary.json`,
`coverage_audit.csv`, and `coverage_audit_summary.json`.
`--bundle-audit` also produces `data/sweap/events_v2/sweap_audit_bundle.zip`,
containing just the compact per-product JSON metadata and CSV audit reports,
not the potentially large NPZ particle arrays.
`--no-resume` forces regeneration. Changed extraction semantics and time
windows automatically invalidate stale resumes.

**Full 599-event run:** only after auditing the one-event results, omit
`--max-events 1` and raise the total download budget to match
available disk and network capacity (e.g. `--max-total-mib 8192`).
`--event-products spc-l2i,spc-l3i` restricts the source products if needed.

**First-event evidence and known source incompleteness (2018-11-01):**
The FIELDS candidate starts at 01:32:05 UTC and ends at 02:00:48 UTC
(28m43s). Its audit bundle showed 13 saved product archives, zero
variables skipped relative to the *source CDFs*, and non-fill event-time
science measurements in all 13. This does **not** imply the original
published product's expected variables are all present. For instance,
the captured SPC L2 time-subset CDF has `Epoch`, current and `mv_lo`
but lacks `measurement_time` and `mv_hi` even though these quantities
are documented in the nominal SPC L2 description; therefore derived
per-voltage-step timestamps are explicitly **unavailable** for this
event. The captured SPA electron 3D subset lacks `THETA`, `PHI`
and `QUALITY_FLAG`; the alpha 3D subset lacks `THETA` and `PHI`.
Those variables must not be reconstructed from guesses. SPAN-B's
electron 3D file does provide angles. Whether these differences reflect
the 2018 native CDF itself or the NASA time-subset transformation needs
a source-versus-subset comparison before making maximum-coverage claims.

The offline `--audit-events` command now inspects saved NPZ quality
flags at exact FIELDS event times, reports their zero/nonzero counts,
and checks for missing science-critical fields. It does not download
new data. The bundle can be re-created after a `git pull` without
re-extracting the event:
`python project.py sweap --audit-events --bundle-audit --event-out data/sweap/events_v2`.

**Scientific caveats:** The 13 selected sources represent eight SWEAP
measurement families, not every possible NASA mission dataset or
unreleased raw Level 0/1 telemetry. L2 is calibrated detector data;
L3 includes derived plasma quantities. Raw spacecraft time tagging,
instrument integration durations, DQF definitions, 2018 viewing geometry,
and channel-dependent systematic uncertainties must still be validated
before inferring physical changes at an individual magnetic switchback.
Some variables can have zero timestamps within the event even with
useful surrounding data. Empty observations are never zero flux.


## Repository guide

```text
project-solar/
├── project.py                    # Stable, ASCII one-command entry point
├── README.md                     # Current overview (this file)
├── The_Research.md               # Scientific objectives and research notes
├── .gitignore                    # Keeps NASA CDFs, caches and outputs local
└── code/
    ├── FIELDS/                   # Magnetic-field data and event detection
    ├── IS⊙IS/                    # exactly isois_reader.py + isois_analyzer.py
    └── SWEAP/                    # Exploratory plasma data inspection
```

The instrument folder is literally `code/IS⊙IS/` on GitHub. Older local
commands using `code/IS_O_IS/` refer to the **previous folder name**.
Prefer `project.py` so older path examples do not break again.

- [Research objectives and scientific background](The_Research.md)

## Data preservation and reproducibility

**`data/` is ignored by Git.** Your 599-event CSV, original NASA CDF files,
compressed per-day particle partitions and coverage report CSV/JSONs remain
on the analysis computer. A `git pull` does not fetch or delete them.
Keep an independent backup of this directory before any restructuring.
Never run `git clean -fdx` or force-reset your working tree to tidy folders.

Run `python project.py test` to execute the eight offline instrument self-tests.
These checks validate software behavior, **not** scientific calibration,
true data availability or causal particle acceleration.

## Scientific interpretation

Missing telemetry or absent particle counts must never be interpreted
as quiet solar wind. Event-selection criteria, actual integration exposure,
detector channels, energy bins, measurement uncertainty and quality flags
must be retained. Hourly measurements provide context; they cannot
reconstruct missing 1-second or 10-second observations.

Primary data come from the Parker Solar Probe FIELDS, IS⊙IS and SWEAP
teams and their official NASA-accessible archives. See
[The_Research.md](The_Research.md) for the study's scientific motivation
and limitations.
