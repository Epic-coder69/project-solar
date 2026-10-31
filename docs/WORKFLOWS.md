# Detailed workflows (preserved technical reference)

**For the current entry point and stable commands, use the [root README](../README.md) and [repository map](REPOSITORY_MAP.md).**
This is the previous long-form README preserved for scientific context and historical commands. Some example commands invoke instrument paths directly; use `python project.py ...` for routine work instead. Status tables in this preserved reference may predate the latest audit. Commands are run from the repository root.

---

# Previous Project Solar workflow reference

**Parker Solar Probe | Magnetic-field structures, energetic particles, and their correlations**

## Research objective

**The central goal of Project Solar is to determine how magnetic-field structures and variability measured by Parker Solar Probe's FIELDS instrument relate to other physical measurements, especially energetic particles measured by IS⊙IS.**

We start with magnetic-field rotations and *switchback candidates* because they provide identifiable time intervals and quantitative magnetic properties. **Detecting switchbacks is a means to the scientific question, not the final result.**

The research asks:

- Do energetic-particle count rates or calibrated intensities change around strong magnetic rotations and radial-field reversals?
- How do particle responses depend on magnetic deflection angle, event duration, radial reversal, and changes in the magnetic-field magnitude?
- Are there energy-, species-, or detector-dependent relationships, and do changes precede or follow the magnetic event?
- Do similar particle variations occur during comparable **non-switchback control intervals**?
- How do solar-wind plasma conditions (SWEAP), electromagnetic variability (FIELDS), heliocentric distance, and data coverage change the interpretation?

Our intended result is a **time-aligned, quality-controlled correlation study** of FIELDS magnetic measurements and complementary PSP observations—not an assumption that every magnetic switchback accelerates particles. Statistical association alone will not establish causation or distinguish local acceleration from transport and other solar-wind processes.

See [The_Research.md](../The_Research.md) for scientific background, research questions, methodology, and previous exploratory work.

## One proposed switchback formation mechanism

![Interchange reconnection: an open magnetic field line interacts with a closed coronal loop, reconnects, and launches an outward-propagating magnetic switchback](https://commons.wikimedia.org/wiki/Special:FilePath/Illustration_of_global_magnetic_field_circulation_enabled_by_interchange_reconnection.jpg)

*Illustration: open and closed coronal magnetic field lines (A) interact, (B) undergo interchange reconnection, and (C) launch an S-shaped magnetic deflection into the corona. This is **one proposed mechanism**, not an established origin of every switchback, and not itself evidence of particle acceleration. Source: [Fisk & Kasper (2020)](https://doi.org/10.3847/2041-8213/ab8acd); image by L. A. Fisk and J. C. Kasper, [Wikimedia Commons](https://commons.wikimedia.org/wiki/File:Illustration_of_global_magnetic_field_circulation_enabled_by_interchange_reconnection.jpg), [CC BY 3.0](https://creativecommons.org/licenses/by/3.0/).*

## Current research status

| Component | Status |
| --- | --- |
| **FIELDS magnetic-field acquisition and rotation screening** | Implemented and run on the November–December 2018 interval |
| **FIELDS detailed magnetic-event characterization** | Implemented; candidate boundaries/physical classification still require independent validation |
| **Strong magnetic-candidate archival** | Completed for the November–December 2018 run, with checksums and a retention manifest |
| **IS⊙IS summary-rate reader** | Initial implementation in GitHub; survey-only, requires real-CDF end-to-end verification |
| **IS⊙IS event-centered rate characterization** | Initial implementation; pre/during/post comparisons against FIELDS event times, not a validated correlation test |
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
           +-- IS⊙IS EPI-Lo / EPI-Hi | 
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
    ├── FIELDS/
    │   ├── magnetic_data_reader.py
    │   ├── magnetic_data_analysis.py
    │   └── external_libraries_requirements.txt
    ├── IS⊙IS/
    │   ├── particle_data_reader.py
    │   └── particle_data_analysis.py
    └── SWEAP/
        └── sweap_data_reader.py
```

- **FIELDS reader:** obtains full-cadence MAG data and screens magnetic investigation windows, with progress reporting and daily checkpoints.
- **FIELDS analysis:** evaluates magnetic rotations, writes event/audit/registry CSVs, selects strong magnetic candidates, and optionally archives verified data before cleaning scoped original-CDF caches.
- **IS⊙IS reader:** reads selected **Level-2 summary count-rate** variables, preserves per-sample integration intervals and units, applies metadata-based invalid-value checks, and writes a CSV and checksum manifest.
- **IS⊙IS analysis:** aligns those survey-rate integrations with existing FIELDS strong-event intervals to summarize exposure-weighted pre/during/post rates and data coverage. This is **exploratory**, not a statistically validated correlation or calibrated flux analysis.

The old experimental work is preserved on the [pre-cleanup archive branch](https://github.com/Epic-coder69/project-solar/tree/archive/pre-cleanup-2026-10-25).

## SWEAP: inspect NASA plasma measurements and actual cadence

The exploratory reader at [`code/SWEAP/sweap_data_reader.py`](../code/SWEAP/sweap_data_reader.py) downloads official **daily CDFs** (cached under Git-ignored `data/sweap/`) or inspects an already-downloaded CDF. It catalogs variable names and units, reports **measured EPOCH spacing** (median/p10/p90/largest gap), shows global quality-flag counts when available, masks fill values and prints the unresampled measurements inside a selectable eight-minute window.

It currently supports `spc-l3i` (**default**: derived proton density/velocity/thermal-speed moments and fits), `spc-l2i` (cup current and charge-flux spectra), `spi-proton-l3` (partial SPAN-I proton moments), and `spe-pad-l3` (electron pitch-angle distribution metadata and cadence). The script is for data exploration—not validated plasma moments, interpretation, or automatic switchback detection.

Run from the repository root:

```bash
python -m pip install -r code/FIELDS/external_libraries_requirements.txt
python code/SWEAP/sweap_data_reader.py --self-test
python code/SWEAP/sweap_data_reader.py --list-products

# Download one NASA day and inspect the first 8 minutes with published records.
python code/SWEAP/sweap_data_reader.py --product spc-l3i \
  --date 2018-11-06 --download

# Replace the example UTC start with an actual FIELDS candidate's start.
python code/SWEAP/sweap_data_reader.py --product spc-l3i \
  --date 2018-11-06 --download --start 2018-11-06T12:00:00 \
  --minutes 8 --csv data/sweap/example_window.csv \
  --json data/sweap/example_window_diagnostics.json

# Inspect an existing NASA CDF without network access.
python code/SWEAP/sweap_data_reader.py --product spc-l3i \
  --file data/sweap/spc-l3i/ACTUAL_FILENAME.cdf --minutes 8
```

**Important:** The nominal SPC archive cadence (~0.2185 s) is not guaranteed real coverage. SWEAP EPOCH records may be irregular and mark the start of an instrument spectrum, not a uniform measurement integration. The tool measures available timestamp separations, not true exposure time; a valid scientific comparison also needs quality flags, uncertainties, field-of-view checks and alignment with the actual FIELDS event interval. SPAN-I moments may be incomplete due to field-of-view limitations. The CDF's native measurement variables must not be confused with raw Level-0 detector packets.

NASA references: [SWEAP data user guide](https://spdf.gsfc.nasa.gov/pub/data/psp/sweap/sweap_data_user_guide.pdf) and [SPC L3 catalog](https://spase-metadata.org/NASA/NumericalData/ParkerSolarProbe/SWEAP/SPC/Level3/SolarWindMomentsFits/PT0.2185S).

## Setup

Run commands from the **project root**. Python 3.10+ recommended:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r code/FIELDS/external_libraries_requirements.txt
```

All NASA CDF downloads, candidate tables, retained archives, and intermediate measurements belong in **`data/`**, which is excluded from Git. Back up that directory separately.

### Self-tests (no NASA downloads)

```bash
python code/FIELDS/magnetic_data_reader.py --self-test
python code/FIELDS/magnetic_data_analysis.py --self-test
python code/IS⊙IS/particle_data_reader.py --self-test
python code/IS⊙IS/particle_data_analysis.py --self-test
```

Passing synthetic tests does **not** prove that the physics, instrument calibrations, real-data coverage, or cross-instrument correlations are correct.

## FIELDS: reproduce the analysis workflow

A short screening example:

```bash
python code/FIELDS/magnetic_data_reader.py --scan \
  --start 2018-11-07T03:30:00 --duration 10m \
  --output data/example_magnetic_candidates.csv
```

The reader saves screened **investigation windows**, not confirmed physical events. For a completed reader scan:

```bash
python code/FIELDS/magnetic_data_analysis.py \
  --windows-csv data/example_magnetic_candidates.csv --no-plots
```

The completed two-month study used:

```bash
python code/FIELDS/magnetic_data_reader.py --scan \
  --start 2018-11-01 --duration 2mo \
  --output data/nov_dec_2018_candidates.csv

python code/FIELDS/magnetic_data_analysis.py \
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

## IS⊙IS: native-cadence particle data for FIELDS events

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

**1-second rates-only coverage:** `SECOND-RATES` contains LET1 A/B, LET2 C,
and HET A/B proton and electron **raw counts and count rates**. It has no
`*_Flux` arrays and is not an energy-resolved calibrated differential-flux
product. It is included by default for precise timing, alongside the 10-second
energy-resolved science products. All measurements retain their detector identity,
original time integration, CDF quality flags, and source metadata.

Helium is not a substitute for the Z=6–28 LET1 heavy-ion channel. The 60 s
species-specific measurements are not identical to the original combined heavy-ion
summary rate. An additional 300 s grouped product is available using
`--products let1-rates300`.

Synthetic tests (no downloads):

```bash
python code/IS⊙IS/particle_data_reader.py --highres --self-test
```

Read the **already downloaded** HET 10-second sample (preserving the raw CDF):

```bash
python code/IS⊙IS/particle_data_reader.py --highres \
  --file data/isois/epihi/psp_isois-epihi_l2-het-rates10_20181106_v26.cdf \
  --events-csv data/full_analysis_nov_dec_2018_candidates_strong_events.csv
```

Plan all four categories against the verified FIELDS switchback candidates:

```bash
python code/IS⊙IS/particle_data_reader.py --highres \
  --events-csv data/full_analysis_nov_dec_2018_candidates_strong_events.csv \
  --plan-only
```

Acquire all candidate-aligned CDFs at their native resolutions (default: 1 s rate-only quicklook, 10 s LET1, 60 s heavy ions and 10 s HET). This may involve many downloads:

```bash
python code/IS⊙IS/particle_data_reader.py --highres \
  --events-csv data/full_analysis_nov_dec_2018_candidates_strong_events.csv \
  --download --resume
```

To validate the 1-second rate product already cached for November 6, 2018:

```bash
python code/IS⊙IS/particle_data_reader.py --highres \
  --start 2018-11-06 --end 2018-11-07 \
  --products second-rates --download --resume
```

### Event-centered observation windows: 10 seconds before/after

High-resolution EPI-Hi extraction now defaults to the exact FIELDS candidate
duration **plus 10 seconds before and 10 seconds after each candidate**.
Do not replace each candidate's duration with a fixed 8-minute interval:
467 seconds is the median of 599 FIELDS candidates, so the median total
observation window is **487 seconds (8 minutes 7 seconds)**.

The reader accepts \`--before-seconds\` and \`--after-seconds\` (both default
to 10). Older explicit \`--before-minutes\` and \`--after-minutes\` arguments
remain supported for exploratory comparisons and override the seconds setting.

Run the **nominal timing coverage audit**, including if downloads are still
incomplete. It scans existing NPZ partitions, keeps detector products distinct,
and verifies stored SHA-256 file manifests by default:

\`\`\`bash
python code/IS⊙IS/particle_coverage_audit.py --self-test
python code/IS⊙IS/particle_coverage_audit.py \
  --events-csv data/full_analysis_nov_dec_2018_candidates_strong_events.csv \
  --before-seconds 10 --after-seconds 10
\`\`\`

The audit writes:
- \`data/isois/highres/coverage_10s.csv\`: each candidate × product, with the
  *nominal unfiltered* coverage fraction before, during and after; complete
  integrations during the event; and two distinct timing-coverage indicators
- \`data/isois/highres/coverage_10s.summary.json\`: product and total counts,
  included source partitions, checksum failures and explicit science caveats

A \`fine_10s_window_temporal_coverage\` result requires a product with at most
10-second median integrations, at least 90% nominal integration coverage in
each event segment, and at least three completely contained integrations
during the event. \`precise_10s_flanks_approx_1s\` imposes tighter sampling:
at most two-second integrations and at least five *complete* integrations in
each 10-second flank. These are **operational timing checks**, not declarations
of scientifically usable exposure, statistical significance or independent
particle counts. In particular, a 60-second integration can overlap a 10-second
flank but cannot isolate it from the switchback.

Until the quality flags, energy-bin validity and Poisson count statistics are
evaluated, **do not claim appropriate scientific coverage merely because the
timing audit passes**. Missing source days are not quiet solar wind.

### Interrupted downloads and safe re-runs

The public IS⊙IS archive may disconnect during a long multi-day download.
The downloader now retries individual CDF transfers, preserves `.part` files
and requests HTTP Range resumes when supported by the server; if Range is
ignored it starts the file again instead of appending duplicates. Each
successfully decoded product/day is saved *immediately*, and
`data/isois/highres/index.json` is checkpointed throughout the run. It
records `status: partial` or `status: complete`, plus missing and failed
product-days. `--resume` validates the NPZ and manifest SHA-256 before reusing
them; it does not require rerunning FIELDS.

```bash
python code/IS⊙IS/particle_data_reader.py --highres \
  --self-test
python code/IS⊙IS/particle_data_reader.py --highres \
  --events-csv data/full_analysis_nov_dec_2018_candidates_strong_events.csv \
  --download --resume --download-attempts 4
```

The code may finish with exit status 2 and an `[INCOMPLETE]` message if any
source product-day is unavailable. This indicates an *explicit missing-data
record*, not that all events were observed or that the downloaded partition
files were lost. A later re-run with `--resume` retries missing or failed
source days.

The reader validates the known `*_Rate` variables and units, stores counts and
count rates separately, and **refuses to reinterpret rates-only data as
differential flux**. Product-specific fill-value masks are stored as
`<variable>__valid_numeric`. These are *not* a substitute for the independent
`Quality_Flag` validation. Coverage CSVs report **unfiltered time overlap**;
they are not claims of usable detector exposure or statistically significant
particle enhancements.

Outputs under `data/isois/highres/` include `index.json`, one
`<source-file>.npz` per product/day, accompanying `.manifest.json`, and
per-file `.coverage.csv` unfiltered time-overlap audits. Event-window NPZ files include
an `_events_<fingerprint>` suffix to prevent full-day reprocessing from
overwriting event-targeted data. `index.json` summarizes the most recent run. Open NPZ using
`np.load(filename, allow_pickle=False)`. Timestamps are UTC Unix
**nanoseconds**. Coverage is time overlap **before interpreting quality flags**.
Do not treat a sample overlapping a switchback as evidence of particle
acceleration; check quality, count statistics, field-of-view, controls, and
possible time delays. Original NASA CDFs and the FIELDS catalog are never deleted.

## IS⊙IS: multi-resolution SAME-channel recovery of missing coverage

This is a **fallback inventory**, not numerical interpolation or a merged
high-cadence flux time series. It ranks **measurements from the same
telescope, A-side and particle species** at their recorded native
integrations, without double-counting the same detections. It never
substitutes HET protons for LET1 protons, helium for heavy ions, or B-side for
A-side; it does not claim that rate and flux energy passbands are equivalent.

**Known 2018 coverage**: our first EPI-Hi 1/10/60-second encounter archive
gave 295/599 switchback candidates with nominal fine timing. The separate
original 3600-second IS⊙IS summary survey had overlapping integrations for
597/599 events, but those are *hourly context*, not event-resolved counts.

The reader now also recognizes native-cadence EPI-Hi
\`let1-rates3600\` and \`het-rates3600\` archives, with optional intermediate
\`let1-rates300\`, \`het-rates60\`, \`het-rates300\` and
\`let2-rates60\`. The intermediate products are only used where their
archived daily files and matching detector channels actually exist.

**Reader compatibility note (2018 hourly HET):** Some official
`het-rates3600` CDFs declare time-dependent auxiliary variables with
**zero actual records**, such as `A_CNO_SECT_Rate` shaped `(0, 2, 25)`
despite the nonempty `Epoch` array. The reader now records these
variables under `skipped_zero_record_variables` in the per-file manifest
instead of rejecting all otherwise-valid hourly data or manufacturing
zero particle counts. Only science variables actually present are eligible
for source matching.

The fallback analyzer explicitly converts timestamps from the old hourly
summary CSV to Unix **nanoseconds**. Without that conversion, pandas versions
that parse timestamps as microseconds may falsely report no hourly coverage.

First run the existing survey CSV through the multiresolution audit **without
downloading anything**:

\`\`\`bash
python code/IS⊙IS/particle_multires_coverage.py --self-test

python code/IS⊙IS/particle_multires_coverage.py \
  --events-csv data/full_analysis_nov_dec_2018_candidates_strong_events.csv \
  --summary-csv data/isois/fields_event_summary_rates.csv
\`\`\`

Acquire the two genuine **hourly EPI-Hi L2 CDF products** as additional
background information for the exact FIELDS event dates. This does **not**
replace the existing high-resolution archives or summary CSV:

\`\`\`bash
python code/IS⊙IS/particle_data_reader.py --highres \
  --events-csv data/full_analysis_nov_dec_2018_candidates_strong_events.csv \
  --products let1-rates3600 het-rates3600 \
  --download --resume --download-attempts 4
\`\`\`

Optionally test additional intermediate products if the official archive
lists them during the relevant dates:

\`\`\`bash
python code/IS⊙IS/particle_data_reader.py --highres \
  --events-csv data/full_analysis_nov_dec_2018_candidates_strong_events.csv \
  --products let1-rates300 het-rates60 het-rates300 let2-rates60 \
  --download --resume --download-attempts 4
\`\`\`

Re-run the multires audit after any download. It creates:

- \`data/isois/highres/multires_best_per_event_channel.csv\`: each candidate
  × each original particle channel, one preferred product/provenance selected
  by *timing tier*, with counts of alternate overlapping source products
- \`data/isois/highres/multires_observed_sources.csv\`: all recorded
  channel-specific overlapping product/event alternatives, never summed
- \`data/isois/highres/multires_source_inventory.csv\`: measured integration
  durations and source-product coverage
- \`data/isois/highres/multires_coverage_summary.json\`: aggregate coverage,
  integrity issues and science caveats

**Coverage-tier interpretation:** \`fine_event_and_flanks\` means nominal
1–10-second timing during the event and in both ±10-second windows.
\`fine_event_only\` means high-cadence event coverage without full flank
coverage. \`coarse_event_scale\` is 60–300-second context with at least
one complete integration inside the event, without resolving the 10-second
flanks. \`hourly_background_only\` is a 3600-second observation, which
cannot establish whether the particle rate changed during that particular
switchback. \`no_observation\` means no matched measurement was found.
All tiers refer to **time availability, not quality-approved physical data**.

The combined heavy-ion \`A_Heavy_Rate_TS\` is a Z=6–28 group. The 60 s
separately identified species and 300 s grouped elemental measurements are
**not the same variable**; the multires preference deliberately does not
claim they are calibrated replacements for the combined group. Their
native-resolution CDFs remain available for a distinct composition analysis.

All EPI-Hi products may have different energy passbands; no energy-bin
remapping, summation or flux/rate conversion is carried out. Numeric
validity masks and quality flags require separate instrument-specific
evaluation, as do Poisson statistics and direction-dependent effects.

## IS⊙IS: initial exploratory survey workflow

First inspect the variables and metadata in an existing **summary** CDF:

```bash
python code/IS⊙IS/particle_data_reader.py \
  --inspect path/to/psp_isois_l2-summary_20181106_v26.cdf
```

For a one-day survey download and extraction:

```bash
python code/IS⊙IS/particle_data_reader.py \
  --download --start 2018-11-06 --end 2018-11-07 \
  --output data/isois/summary_rates_20181106.csv
```

Then run the **event-centered survey comparison** against the previously retained FIELDS strong-event catalog:

```bash
python code/IS⊙IS/particle_data_analysis.py \
  --rates-csv data/isois/summary_rates_20181106.csv \
  --events-csv data/full_analysis_nov_dec_2018_candidates_strong_events.csv \
  --output data/isois/fields_event_particle_survey.csv
```

With only one day of particle data, **most FIELDS events will correctly have insufficient particle coverage**. The analysis writes per-event/per-rate-channel pre/during/post estimates and explicitly flags inadequate exposure; it does not claim a statistical correlation.

**Crucial limitation:** Level-2 IS⊙IS **summary count rates are not calibrated differential particle fluxes**. This starter workflow is suitable for exploring candidate time intervals only. Science-grade particle intensity/spectral comparisons require the appropriate EPI-Lo and EPI-Hi science products, energy-bin and species metadata, proper uncertainties and quality flags, plus independently chosen quiet/control intervals.

## Scientific interpretation and next steps

1. **Validate** the FIELDS candidate selection and the IS⊙IS reader against real observations and instrument documentation.
2. **Acquire science-grade particle products** and calculate physically interpretable energy- and species-resolved measurements.
3. **Build the correlation layer**: match exposure intervals, test pre/during/post differences and possible lags, compare with non-switchback controls, account for unequal coverage and changing spacecraft/solar-wind conditions, report uncertainty and significance.
4. **Extend the measurements** to SWEAP plasma parameters and additional FIELDS electric/wave diagnostics where data quality permits.

Even a reproducible temporal coincidence does not by itself prove **particle acceleration inside a switchback**; particle transport and common solar-wind conditions must also be considered.

## Data and instrument acknowledgments

- [PSP/FIELDS data portal](https://fields.ssl.berkeley.edu/data/) and NASA CDAWeb.
- Parker Solar Probe IS⊙IS instrument team and public Level-2 products.
- Parker Solar Probe SWEAP team for any future plasma or ephemeris products.
- [Martinović et al. switchback catalog](https://doi.org/10.5281/zenodo.4579821) for exploratory external comparison.

**Project Solar is a work in progress. The scientific objective is correlation and physical interpretation; automatic event detection is one necessary supporting step.**
