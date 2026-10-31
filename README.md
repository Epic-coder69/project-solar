# Project Solar

**Parker Solar Probe | Magnetic-field structures, energetic particles, and their correlations**

## Research objective

**The central goal of Project Solar is to determine how magnetic-field structures and variability measured by Parker Solar Probe's FIELDS instrument relate to other physical measurements, especially energetic particles measured by ISʘIS.**

We start with magnetic-field rotations and *switchback candidates* because they provide identifiable time intervals and quantitative magnetic properties. **Detecting switchbacks is a means to the scientific question, not the final result.**

The research asks:

- Do energetic-particle count rates or calibrated intensities change around strong magnetic rotations and radial-field reversals?
- How do particle responses depend on magnetic deflection angle, event duration, radial reversal, and changes in the magnetic-field magnitude?
- Are there energy-, species-, or detector-dependent relationships, and do changes precede or follow the magnetic event?
- Do similar particle variations occur during comparable **non-switchback control intervals**?
- How do solar-wind plasma conditions (SWEAP), electromagnetic variability (FIELDS), heliocentric distance, and data coverage change the interpretation?

Our intended result is a **time-aligned, quality-controlled correlation study** of FIELDS magnetic measurements and complementary PSP observations—not an assumption that every magnetic switchback accelerates particles. Statistical association alone will not establish causation or distinguish local acceleration from transport and other solar-wind processes.

See [The_Research.md](The_Research.md) for scientific background, research questions, methodology, and previous exploratory work.

## One proposed switchback formation mechanism

![Interchange reconnection: an open magnetic field line interacts with a closed coronal loop, reconnects, and launches an outward-propagating magnetic switchback](https://commons.wikimedia.org/wiki/Special:FilePath/Illustration_of_global_magnetic_field_circulation_enabled_by_interchange_reconnection.jpg)

*Illustration: open and closed coronal magnetic field lines (A) interact, (B) undergo interchange reconnection, and (C) launch an S-shaped magnetic deflection into the corona. This is **one proposed mechanism**, not an established origin of every switchback, and not itself evidence of particle acceleration. Source: [Fisk & Kasper (2020)](https://doi.org/10.3847/2041-8213/ab8acd); image by L. A. Fisk and J. C. Kasper, [Wikimedia Commons](https://commons.wikimedia.org/wiki/File:Illustration_of_global_magnetic_field_circulation_enabled_by_interchange_reconnection.jpg), [CC BY 3.0](https://creativecommons.org/licenses/by/3.0/).*

## Current research status

| Component | Status |
| --- | --- |
| **FIELDS magnetic-field acquisition and rotation screening** | Implemented and run on the November–December 2018 interval |
| **FIELDS detailed magnetic-event characterization** | Implemented; candidate boundaries/physical classification still require independent validation |
| **Strong magnetic-candidate archival** | Completed for the November–December 2018 run, with checksums and a retention manifest |
| **ISʘIS summary-rate reader** | Initial implementation in GitHub; survey-only, requires real-CDF end-to-end verification |
| **ISʘIS event-centered rate characterization** | Initial implementation; pre/during/post comparisons against FIELDS event times, not a validated correlation test |
| **Calibrated EPI-Lo / EPI-Hi flux and energy-spectral analysis** | Not yet implemented |
| **Matched controls, statistical tests, and final cross-instrument correlations** | Not yet implemented |
| **SWEAP plasma and additional FIELDS wave/electric-field correlations** | Planned extensions |

### FIELDS survey: November 1–December 31, 2018

The completed **61-calendar-day requested interval** produced:

| Metric | Recorded result |
| --- | ---: |
| Reader investigation windows | **866** |
| Preliminary magnetic rotations characterized | **3,701** |
| Magnetically selected strong candidates | **599** |
| Valid 1-second magnetic-vector coverage of the requested UTC interval | **68.4%** |

These are **magnetic candidates, not 599 physically confirmed switchbacks**. Missing telemetry is not quiet solar wind; the observation period is not 61 uninterrupted days of measurements.

The strong-candidate screen requires a native-cadence peak deflection of at least **90°**, a radial-field reversal, a median absolute magnetic-magnitude deviation of at most **10%** from the chosen reference, uncensored event boundaries, and valid native measurements. These are **operational selection criteria**, not independent physical verification.

The associated data products were generated on the research computer under the Git-ignored `data/` directory. They are **not stored in this GitHub repository**.

## How the research pipeline fits together

```text
Parker Solar Probe observations
           |
           +-- FIELDS MAG, RTN magnetic vectors
           |       |
           |       +-- magnetic_data_reader.py
           |       |      Download, verify, screen time intervals
           |       |
           |       +-- magnetic_data_analysis.py
           |              Detect/characterize rotations, assign event IDs,
           |              audit data quality, retain strong magnetic events
           |                          |
           |                     FIELDS event catalog
           |                          |
           +-- ISʘIS EPI-Lo / EPI-Hi | 
           |       |                 |
           |       +-- particle_data_reader.py
           |       |      Current: summary-rate survey extraction
           |       |
           |       +-- particle_data_analysis.py
           |              Current: exposure-aware pre/during/post
           |              survey rates at FIELDS event times
           |                          |
           |                     Particle measurements
           |                          |
           +--------------------------+
                                      |
                              CORRELATION ANALYSIS
                           (planned; not implemented)
                                      |
                    Matched timing, controls, uncertainties,
                    time lags, energy/species dependence,
                    and scientific interpretation
```

**Common data requirements:** preserve UTC timing, each instrument's actual integration interval/cadence, detector and channel metadata, units, valid-data masks, coverage, provenance, and stable magnetic `event_id` values. A missing observation is not a zero particle rate.

## Repository structure

```text
Project Solar/
├── README.md
├── The_Research.md
└── code/
    ├── M_FIELD/
    │   ├── magnetic_data_reader.py
    │   ├── magnetic_data_analysis.py
    │   └── external_libraries_requirements.txt
    └── IS_O_IS/
        ├── particle_data_reader.py
        └── particle_data_analysis.py
```

- **FIELDS reader:** obtains full-cadence MAG data and screens magnetic investigation windows, with progress reporting and daily checkpoints.
- **FIELDS analysis:** evaluates magnetic rotations, writes event/audit/registry CSVs, selects strong magnetic candidates, and optionally archives verified data before cleaning scoped original-CDF caches.
- **ISʘIS reader:** reads selected **Level-2 summary count-rate** variables, preserves per-sample integration intervals and units, applies metadata-based invalid-value checks, and writes a CSV and checksum manifest.
- **ISʘIS analysis:** aligns those survey-rate integrations with existing FIELDS strong-event intervals to summarize exposure-weighted pre/during/post rates and data coverage. This is **exploratory**, not a statistically validated correlation or calibrated flux analysis.

The old experimental work is preserved on the [pre-cleanup archive branch](https://github.com/Epic-coder69/project-solar/tree/archive/pre-cleanup-2026-10-25).

## Setup

Run commands from the **project root**. Python 3.10+ recommended:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r code/M_FIELD/external_libraries_requirements.txt
```

All NASA CDF downloads, candidate tables, retained archives, and intermediate measurements belong in **`data/`**, which is excluded from Git. Back up that directory separately.

### Self-tests (no NASA downloads)

```bash
python code/M_FIELD/magnetic_data_reader.py --self-test
python code/M_FIELD/magnetic_data_analysis.py --self-test
python code/IS_O_IS/particle_data_reader.py --self-test
python code/IS_O_IS/particle_data_analysis.py --self-test
```

Passing synthetic tests does **not** prove that the physics, instrument calibrations, real-data coverage, or cross-instrument correlations are correct.

## FIELDS: reproduce the analysis workflow

A short screening example:

```bash
python code/M_FIELD/magnetic_data_reader.py --scan \
  --start 2018-11-07T03:30:00 --duration 10m \
  --output data/example_magnetic_candidates.csv
```

The reader saves screened **investigation windows**, not confirmed physical events. For a completed reader scan:

```bash
python code/M_FIELD/magnetic_data_analysis.py \
  --windows-csv data/example_magnetic_candidates.csv --no-plots
```

The completed two-month study used:

```bash
python code/M_FIELD/magnetic_data_reader.py --scan \
  --start 2018-11-01 --duration 2mo \
  --output data/nov_dec_2018_candidates.csv

python code/M_FIELD/magnetic_data_analysis.py \
  --windows-csv data/nov_dec_2018_candidates.csv \
  --no-plots --keep-strong-only
```

**Do not rerun the two-month scan simply to inspect the existing results.** It may download large volumes of NASA data. The already-produced catalogs should be used for subsequent particle comparisons.

Key local outputs are:

```text
data/nov_dec_2018_candidates.csv
data/full_analysis_nov_dec_2018_candidates_events.csv
data/full_analysis_nov_dec_2018_candidates_audit.csv
data/full_analysis_nov_dec_2018_candidates_registry.csv
data/full_analysis_nov_dec_2018_candidates_strong_events.csv
data/full_analysis_nov_dec_2018_candidates_strong_registry.csv
data/full_analysis_nov_dec_2018_candidates_strong_archive_index.csv
data/full_analysis_nov_dec_2018_candidates_strong_retention.json
data/full_analysis_nov_dec_2018_candidates_strong_archive/
```

The high-resolution archives retain quality-screened magnetic vectors and context; they are **not bit-for-bit replacements for NASA's original CDFs and full instrument metadata**. The opt-in `--keep-strong-only` option can delete eligible local raw CDFs **only after** completed-scan, audit, and archive verification; missing NASA telemetry alone does not veto a verified run. Use `--finalize-existing` with the other required flags to complete archival from an already-finished analysis without repeating computation.

## ISʘIS: native-cadence particle data for FIELDS events

The **high-resolution** reader is separate from the hourly summary-survey mode.
It retains native time integrations, energy bins, counts, differential flux,
uncertainties, pointing/quality data, and all supported science-array dimensions
in per-day compressed NPZ files. **It does not resample the particle data.**

| Original summary channel | Higher-resolution source | Native cadence |
| --- | --- | --- |
| A_H_Rate_TS (LET1 protons) | LET1-RATES10 | 10 s |
| A_Heavy_Rate_TS (Z 6–28 LET1 heavy-ion group) | LET1-RATES60 individual heavy-ion species (C, N, O, Ne, Si, Fe, etc.) | 60 s; **not the same combined quantity** |
| HET_A_H_Rate_TS | HET-RATES10 | 10 s |
| HET_A_Electrons_Rate_TS | HET-RATES10 | 10 s |

`SECOND-RATES` is included for **1-second proton/electron quicklook**, including with the default high-resolution product selection,
but it has reduced energy-bin detail and is **not a substitute for the 10-second
science products for quantitative studies**. Neither He nor separately identified heavy elements are interchangeable with the original Z=6–28 summary count-rate channel. The 300 s grouped CNO/Ne–Si/Fe directional product is also supported via `--products let1-rates300`, but is not the default.

Synthetic tests (no downloads):

```bash
python code/IS_O_IS/particle_data_reader.py --highres --self-test
```

Read the **already downloaded** HET 10-second sample (preserving the raw CDF):

```bash
python code/IS_O_IS/particle_data_reader.py --highres \
  --file data/isois/epihi/psp_isois-epihi_l2-het-rates10_20181106_v26.cdf \
  --events-csv data/full_analysis_nov_dec_2018_candidates_strong_events.csv
```

Plan all four categories against the verified FIELDS switchback candidates:

```bash
python code/IS_O_IS/particle_data_reader.py --highres \
  --events-csv data/full_analysis_nov_dec_2018_candidates_strong_events.csv \
  --plan-only
```

Run the native-cadence archive acquisition only **after** real CDF validation (default: 10 s LET1, 60 s LET1 heavy ions, 10 s HET, and optional-for-physics 1 s quicklook):

```bash
python code/IS_O_IS/particle_data_reader.py --highres \
  --events-csv data/full_analysis_nov_dec_2018_candidates_strong_events.csv \
  --download --resume
```

Outputs under `data/isois/highres/` include `index.json`, one
`<source-file>.npz` per product/day, accompanying `.manifest.json`, and
per-file `.coverage.csv` time-overlap audits. Open NPZ using
`np.load(filename, allow_pickle=False)`. Timestamps are UTC Unix
**nanoseconds**. Coverage is time overlap **before interpreting quality flags**.
Do not treat a sample overlapping a switchback as evidence of particle
acceleration; check quality, count statistics, field-of-view, controls, and
possible time delays. Original NASA CDFs and the FIELDS catalog are never deleted.

## ISʘIS: initial exploratory survey workflow

First inspect the variables and metadata in an existing **summary** CDF:

```bash
python code/IS_O_IS/particle_data_reader.py \
  --inspect path/to/psp_isois_l2-summary_20181106_v26.cdf
```

For a one-day survey download and extraction:

```bash
python code/IS_O_IS/particle_data_reader.py \
  --download --start 2018-11-06 --end 2018-11-07 \
  --output data/isois/summary_rates_20181106.csv
```

Then run the **event-centered survey comparison** against the previously retained FIELDS strong-event catalog:

```bash
python code/IS_O_IS/particle_data_analysis.py \
  --rates-csv data/isois/summary_rates_20181106.csv \
  --events-csv data/full_analysis_nov_dec_2018_candidates_strong_events.csv \
  --output data/isois/fields_event_particle_survey.csv
```

With only one day of particle data, **most FIELDS events will correctly have insufficient particle coverage**. The analysis writes per-event/per-rate-channel pre/during/post estimates and explicitly flags inadequate exposure; it does not claim a statistical correlation.

**Crucial limitation:** Level-2 ISʘIS **summary count rates are not calibrated differential particle fluxes**. This starter workflow is suitable for exploring candidate time intervals only. Science-grade particle intensity/spectral comparisons require the appropriate EPI-Lo and EPI-Hi science products, energy-bin and species metadata, proper uncertainties and quality flags, plus independently chosen quiet/control intervals.

## Scientific interpretation and next steps

1. **Validate** the FIELDS candidate selection and the ISʘIS reader against real observations and instrument documentation.
2. **Acquire science-grade particle products** and calculate physically interpretable energy- and species-resolved measurements.
3. **Build the correlation layer**: match exposure intervals, test pre/during/post differences and possible lags, compare with non-switchback controls, account for unequal coverage and changing spacecraft/solar-wind conditions, report uncertainty and significance.
4. **Extend the measurements** to SWEAP plasma parameters and additional FIELDS electric/wave diagnostics where data quality permits.

Even a reproducible temporal coincidence does not by itself prove **particle acceleration inside a switchback**; particle transport and common solar-wind conditions must also be considered.

## Data and instrument acknowledgments

- [PSP/FIELDS data portal](https://fields.ssl.berkeley.edu/data/) and NASA CDAWeb.
- Parker Solar Probe ISʘIS instrument team and public Level-2 products.
- Parker Solar Probe SWEAP team for any future plasma or ephemeris products.
- [Martinović et al. switchback catalog](https://doi.org/10.5281/zenodo.4579821) for exploratory external comparison.

**Project Solar is a work in progress. The scientific objective is correlation and physical interpretation; automatic event detection is one necessary supporting step.**
