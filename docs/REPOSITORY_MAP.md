# Repository map and safe workflow

This document maps the **actual committed tree**, the role of each script,
and the local-only data directories. It intentionally makes no destructive
changes to the research pipeline.

## Source of truth: what GitHub contains

```text
project-solar/
├── README.md                         # Entry point and detailed workflows
├── The_Research.md                   # Scientific questions and methodology
├── .gitignore                        # Keeps large data, caches, venv out of Git
└── code/
    ├── FIELDS/
    │   ├── magnetic_data_reader.py   # Acquire/select MAG RTN data
    │   ├── magnetic_data_analysis.py # Magnetic candidate detection and audit
    │   └── external_libraries_requirements.txt
    ├── IS⊙IS/
    │   ├── particle_data_reader.py   # Hourly survey reader and --highres CLI
    │   ├── particle_data_analysis.py # Original exploratory survey summaries
    │   ├── particle_highres_reader.py # Native-cadence EPI-Hi CDF/NPZ pipeline
    │   ├── particle_coverage_audit.py # 1/10/60 s timing-coverage check
    │   ├── particle_multires_coverage.py # Same-channel fallback tier audit
    │   └── .gitkeep
    └── SWEAP/
        └── sweap_data_reader.py      # Exploratory plasma-product inspector
```

**Path warning:** the GitHub directory is literally `code/IS⊙IS/` (Unicode
sun-dot). Earlier local terminal examples sometimes used `code/IS_O_IS/`,
which is *not the same pathname*. Do not assume they are interchangeable or
copy one over the other. Confirm the folder on your own computer first:

```bash
cd ~/Desktop/project-solar
git status --short
find code -maxdepth 2 -name 'particle_data_reader.py' -print
git ls-files 'code/*/particle_data_reader.py'
```

For commands below replace `code/IS⊙IS` with your verified local folder
only when necessary. Do **not** rename either folder while an extraction is
running. If both spellings exist locally, inspect `git status` and file
differences before deleting or consolidating anything.

## Data is local, not lost

The Git-ignored `data/` directory stores the FIELDS strong candidate catalog,
its retention/checksum manifest, original NASA CDF files, cached EPI-Hi products,
per-day NPZ partitions and coverage reports. A `git pull` updates **code
and tracked documentation**, not this scientific data. VS Code should display
the local `data/` folder even if it is ignored by Git.

Important local-only inputs and outputs:

| Path | Purpose |
| --- | --- |
| `data/full_analysis_nov_dec_2018_candidates_strong_events.csv` | 599 FIELDS strong magnetic candidates; event boundaries |
| `data/full_analysis_nov_dec_2018_candidates_strong_retention.json` | Catalog retention / integrity manifest |
| `data/isois/fields_event_summary_rates.csv` | Old hourly summary count-rate survey |
| `data/isois/highres/` | EPI-Hi native-cadence NPZ partitions and provenance |
| `data/isois/highres/coverage_10s.csv` | Fine-cadence nominal time-coverage rows |
| `data/isois/highres/multires_best_per_event_channel.csv` | One preferred source per event and original channel |
| `data/isois/highres/multires_observed_sources.csv` | Alternative overlapping sources, not summed |
| `data/isois/highres/multires_coverage_summary.json` | Cross-cadence timing-tier counts and caveats |

A clean-looking GitHub file tree **does not imply the raw data are backed up**.
Back up `data/` separately, including FIELDS strong-event manifests and
downloaded CDF/NPZ files.

## Pipeline roles — what to run when

1. **FIELDS** defines and archives magnetic candidate start/end times.
   Treat these as *candidates*, not independently confirmed physical events.
2. **IS⊙IS acquisition** preserves original particle measurements and their
   actual time integrations. The hourly survey and native EPI-Hi readers are
   intentionally different, not alternative files to arbitrarily combine.
3. **Coverage audits** test *whether and when* measurements exist. The fine
   audit measures 10-second before/after timing; the multires audit reports
   detector/species-specific alternatives and distinguishes hourly background
   from event-resolution data.
4. **Science validation (still pending)** must examine quality flags,
   exposure, energy-bin comparability, count statistics and comparison/control
   event selection before drawing physical conclusions.
5. **SWEAP** presently offers exploratory plasma data inspection, not a
   completed switchback correlation engine.

### Non-destructive check of the current analysis

Activate your environment and use the real local instrument folder:

```bash
cd ~/Desktop/project-solar
source .venv/bin/activate
find code -maxdepth 2 -name 'particle_multires_coverage.py' -print
# Example for GitHub's canonical folder:
python 'code/IS⊙IS/particle_multires_coverage.py' --self-test
python 'code/IS⊙IS/particle_multires_coverage.py' \
  --events-csv data/full_analysis_nov_dec_2018_candidates_strong_events.csv \
  --summary-csv data/isois/fields_event_summary_rates.csv
```

These commands audit existing partitions and do not delete raw CDFs. Do not
rerun the two-month FIELDS scan merely to regenerate particle coverage.

## Latest verified coverage interpretation

The 2018 November–December user-run multires audit reported:

- **599** FIELDS magnetic candidates in the strong-event catalog.
- **295** with nominal 1–10 s IS⊙IS event + 10 s flank timing coverage.
- **302** additional candidates with **hourly particle context only** for
  HET A protons/electrons; hourly averages do not resolve those events.
- **2** without any of the audited original channels.
- LET1 A protons: 295 fine; 301 hourly-only; 1 partial-coarse; 2 missing.
- LET1 combined Z=6–28 heavy-ion rate: 597 hourly context; **no
  event-resolved same-channel equivalent has been established**.

These are **nominal time-coverage figures, not quality-approved particle
detections**. Some high-energy particle counts are sparse, and quality flags,
energy ranges and uncertainty propagation still require independent review.

## What a future structural refactor should do (not yet applied)

- Choose **one portable ASCII package path** (e.g. `code/isois/`) after
  reconciling the Unicode GitHub path with any local aliases.
- Move reusable timing/CDF/manifest utilities into shared modules only after
  tests exist; maintain stable CLI wrappers during migration.
- Move inline self-tests into `tests/` and run them on CI.
- Keep the root README short; place per-instrument details under `docs/`.
- Version lightweight, reproducible *summaries* with explicit provenance when
  appropriate; keep CDFs, huge NPZ arrays and private working data out of Git.

**Safety rule:** do not `git clean -fdx`, delete `data/`, or force-reset a
working tree to fix presentation problems. Inspect `git status` and make a
separate backup of local measurements before any folder rename.
