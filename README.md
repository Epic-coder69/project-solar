# Project Solar

**Parker Solar Probe — magnetic fields, candidate switchbacks, and energetic particles**

This repository contains the **working FIELDS analysis** for Project Solar. Its
magnetic switchback detector is a **candidate-screening method**, not a validated
physical event classifier. Event segmentation is still under investigation.

![Magnetic switchback illustration](figures/switchback_illustration.svg)

## Repository

| Path | Purpose |
|---|---|
| [The_Research.md](The_Research.md) | Research question, methodology and results |
| `code/external_libraries_requirements.txt` | Required Python packages |
| `code/magnetic_data_reader.py` | Read full-cadence FIELDS data and pre-screen arbitrary UTC time intervals; retains lower-cadence readers |
| `code/magnetic_data_analysis.py` | Candidate switchback detector; local or exploratory SPC reference |
| `results/` | Small set of retained detector catalogs and final validation summary |
| `figures/` | Illustrations and figures used in the research document |

**Experimental scripts, raw comparisons, intermediate results and review images
were removed from `main` to reduce clutter. Nothing was destroyed:** the complete
pre-cleanup snapshot remains on
[`archive/pre-cleanup-2026-10-25`](https://github.com/Epic-coder69/project-solar/tree/archive/pre-cleanup-2026-10-25).

## Setup

Run the following commands from the repository root. Python 3.10+ recommended.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r code/external_libraries_requirements.txt
```

Raw NASA CDF data belong in the git-ignored `data/` cache. The loaders use
NASA CDAWeb and fetch missing files when needed; FIELDS high-resolution data
can be large.

## Legacy FIELDS data access

```bash
python code/magnetic_data_reader.py --start 2018-11-06 --days 1 --resolution 1min
python code/magnetic_data_reader.py --start 2018-11-06 --days 1 --resolution 4cyc
```

Legacy products are available using `get_day(day, resolution="1min")`
or `get_day(day, resolution="4cyc")`. For the current full-cadence
workflow, the analysis script generates magnetic event plots directly;
there is no separate visualization script.

## Scan any time range at full cadence

The reader can pre-screen **one continuous UTC time interval**: minutes,
hours, days, weeks, or months. You do not have to choose the six-hour
CDF files or issue separate commands for each date. NASA publishes
full-cadence MAG RTN measurements in six-hour files; the reader retrieves
only the overlapping segments plus the 100-minute context needed to
estimate the magnetic background. Full-cadence CDFs are cached locally
under ignored data/fields_full/ and never discarded when nothing is detected.

The internal work is streaming and memory-bounded. A 1-second representation
is used for *permissive screening*, preserving original full-resolution
samples on disk for detailed analysis. Candidate windows include five
minutes of surrounding context. Missing or sparse coverage is reported;
these are **not physically confirmed switchback events**.

Run synthetic tests without NASA downloads:

```bash
python code/magnetic_data_reader.py --self-test
```

Request any UTC interval in a single command:

```bash
# Ten minutes
python code/magnetic_data_reader.py --scan --start 2018-11-07T03:30:00 --duration 10m

# Two weeks
python code/magnetic_data_reader.py --scan --start 2018-11-01 --duration 2w

# Three calendar months
python code/magnetic_data_reader.py --scan --start 2018-11-01 --duration 3mo

# Exact start and exclusive end
python code/magnetic_data_reader.py --scan --start 2018-11-06T12:00:00 --end 2018-11-06T15:00:00
```

Output goes into a **single candidate CSV** under data/ for the
requested interval; use --output to set another location. **The reader
updates that CSV after every completed day**, while the following days are
still downloading. A small `.checkpoint.json` sidecar keeps unmerged
candidate windows for resuming an interrupted scan:

```bash
# Start a month-long scan; inspect its CSV while the program continues
python code/magnetic_data_reader.py --scan --start 2018-11-01 --duration 1mo

# If interrupted by Ctrl+C, continue from the last completed day
python code/magnetic_data_reader.py --scan --start 2018-11-01 --duration 1mo --resume
```

The same start/end range and --output path must be used when resuming.
The checkpoint sidecar is removed automatically on successful completion.
Partial days currently in progress are **not** checkpointed; the reader
reprocesses that day if interrupted. Original NASA CDF downloads remain cached.
Scans made before this update cannot be resumed retroactively.

The output contains *preliminary investigation windows*, not scientifically
confirmed switchbacks, and the existing event-window merging can still
combine separate magnetic rotations. There is no fixed duration limit,
but missing NASA data, download sizes, disk space and computing time apply.

**Current division of work:** the reader screens original full-cadence
measurements for candidate observation windows. The analysis script now
consumes those candidate CSVs, uses one-second vectors derived from the
full-cadence source to delineate individual magnetic rotations, and revisits
the native-rate observations to calculate physical magnetic metrics. These
events are **not confirmed physical switchbacks**. Independent plasma
verification and synchronized particle visualization remain future work.


## Keep only the strongest events and clear the raw cache

No third script is needed. Run the **reader for the entire selected interval**,
then run the **analysis** with `--keep-strong-only`. This option first writes
the normal event/registry/audit CSVs and then:

1. Selects events with **native peak rotation >=90 degrees**, **radial
   reversal**, **median absolute |B| deviation <=10%** from the reference,
   **uncensored boundaries**, and at least one valid native sample. This is
   the same conservative operational screen used for the 21 promising events
   found in our November 2, 2018 review.
2. Saves `*_strong_events.csv` and `*_strong_registry.csv`, as well as a
   checksummed `*_strong_archive_index.csv` and retention JSON manifest.
3. Preserves **only selected events' quality-screened original-cadence RTN
   magnetic measurements** (event +/- 5 minutes) as gzip CSV files. It also
   keeps their 1-second waveforms, local 181-minute background vectors,
   scalar reference, and deflection angle, for later plotting.
4. Verifies each gzip archive and checksum **before deletion**. Raw deletion
   is skipped if the reader checkpoint is missing/incomplete, scan coverage
   is below 90%, any candidate window failed analysis, or a needed
   background day is missing.
5. Once everything passes, deletes only canonical six-hour original
   `psp_fld_l2_mag_rtn_*_v*.cdf` files directly in
   `data/fields_full/`, wholly within the reader's requested scan range.
   It does not touch other FIELDS resolutions, SWEAP, IS☉IS, source code,
   event catalogs, or files outside the matching scan.

```bash
# First verify synthetic tests (no NASA downloads)
python code/magnetic_data_analysis.py --self-test

# Scan two months (2018-11-01 to 2019-01-01 UTC)
python code/magnetic_data_reader.py --scan --start 2018-11-01 \
  --duration 2mo --output data/nov_dec_2018_candidates.csv

# Analyze all reader candidates; KEEP ONLY strong high-cadence archives;
# automatically clean eligible original CDFs after the full run succeeds.
python code/magnetic_data_analysis.py \
  --windows-csv data/nov_dec_2018_candidates.csv \
  --no-plots --keep-strong-only
```

**Important:** Cleanup happens **at the end of analysis**, not while the
reader runs. Therefore the reader still needs enough disk space to cache up to
two months of six-hour raw files before deletion begins. The analysis does
not resume mid-run; partial output CSVs persist, but interrupted runs need
to be restarted. The raw CDFs are re-downloadable from NASA, whereas event
gzip archives are derived, quality-screened vectors, **not a lossless copy**
of the entire original CDF (instrument metadata and original quality-flag
history are not retained). The small all-event CSVs remain as an audit trail,
but only **strong events** receive high-resolution archives. The option is
opt-in because it actually deletes local cache data.

## Magnetic-event analysis

**Default: analyze the original full-cadence source using the reader's CSV.**

```bash
# Synthetic tests for both full-cadence and legacy logic (no downloads)
python code/magnetic_data_analysis.py --self-test

# First try one targeted candidate (CSV row 11), without plotting
python code/magnetic_data_analysis.py --windows-csv data/november_02_candidates.csv --candidate-index 11 --no-plots

# Optional: retrieve SPC position files when they're not already cached
python code/magnetic_data_analysis.py --windows-csv data/november_02_candidates.csv --candidate-index 11 --no-plots --position-source spc

# Analyze all 20 Nov 2 investigation windows, with figures
python code/magnetic_data_analysis.py --windows-csv data/november_02_candidates.csv
```

The full pipeline uses the reader's `start,end,category` columns, plus
`core_start,core_end` when present. The latter identify which detected
sub-events belong to the screening window, **not** the event boundaries.
The analysis independently finds 45°/30° magnetic rotations lasting at
least 10 seconds, flags incomplete edges, and removes duplicate detections.
It creates:

- `data/full_analysis_november_02_candidates_events.csv` — per-event time boundaries, angle
  from the 1-second derived waveform, native-cadence peak deflection,
  radial reversal, native sample counts, and deviations of `|B|` from
  the local scalar baseline.
- `data/full_analysis_november_02_candidates_audit.csv` — how many events each screening window
  yielded, valid-data/reference coverage, missing background days and errors.
- `data/full_analysis_november_02_candidates_registry.csv` — stable event
  IDs (for unchanged boundaries), event start/peak/end/midpoint UTC times,
  before/during/after windows for instrument correlations, quality and
  magnetic provenance, and spacecraft position when available.
- `figures/full_analysis_*.png` — **1-second display plots** with magnetic
  components, scalar magnitude, local reference and threshold crossings.

The events and registry CSVs share an `event_id` for future joins with
SWEAP plasma, IS☉IS particles and FIELDS wave data. Pre/during/post
windows use five minutes before/after the event, plus a wider 30-minute
context. Correlation intervals are half-open in UTC: the 1-second event
end is treated as inclusive, with an additional `end_exclusive` timestamp
for clean alignment. The registry does not fetch plasma density,
temperature, or particle flux; future modules will use these timestamps
to align instrument data at their actual cadences.

Spacecraft metadata is sourced **only** from `sc_pos_HCI` in
`PSP_SWP_SPC_L3I`; it has HCI Cartesian coordinates in km and derived
distance/latitude/longitude in that same inertial frame (**not** Carrington
coordinates). By default `--position-source cached-spc` uses only locally
cached SPC CDFs and never downloads more data. Specify `--position-source spc`
to obtain missing source CDFs, or `--position-source none` to disable it.
Positions are matched to the event midpoint only when a sample exists within
120 seconds. Otherwise fields remain empty and `position_status` records
why (no fabricated interpolation). Physical plasma-variable quality flags
are not used to mask the spacecraft ephemeris.

These CSVs are updated **after every processed candidate window**, not
only at the end.
The filenames derive from the input CSV and optional `--top` or
`--candidate-index` selection, preventing
scans of different days from overwriting each other's results. The native-cadence magnitude-deviation median uses a
streaming histogram (approximately 0.1 percentage-point precision) to
avoid storing millions of samples in memory. Event timestamps currently
have 1-second precision; the higher-cadence pass measures properties but
does not independently refine sub-second boundaries. Cached NASA CDFs
are reused, but processing large intervals can still be substantial.

**Legacy benchmark:** the 20 frozen survey windows and older 1-minute /
4-samples-per-cycle analysis remain available explicitly:

```bash
python code/magnetic_data_analysis.py --mode legacy --top 1 --reference local --no-plots
python code/magnetic_data_analysis.py --mode legacy --windows-csv data/survey_selected.csv --top 0
```

The experimental `--reference parker-spc` option is **legacy only**.
It uses SWEAP SPC and is not identical to the published SPAN-i-based method.
Legacy results retain their previous `corrected_*` paths.

### Scientific status

- Reading and preprocessing FIELDS magnetic-field vectors: **working**.
- Automatic rotation/switchback *candidate* extraction: **working**.
- New full-cadence analysis with magnetic diagnostics and a correlation-ready
  event registry: **implemented; synthetic self-tests included, requires
  real-data validation on the research Mac**.
- Independent validation: **completed exploratory comparison, inconclusive
  detector ranking**; see [holdout findings](results/holdout_findings.md).
- Physically confirmed switchbacks and robust event boundaries: **not yet**.
- Energetic-particle correlation analysis with IS☉IS: **future work**.

**Reproducibility notice:** The high-resolution loading function has been merged into `code/magnetic_data_reader.py`.
The original local `switchback_survey.py` was not tracked in Git; thus
reproduction of the original event catalogs still needs checking against
actual NASA CDF data on the research Mac. The refactor was tested with
simulated CDF arrays but not yet rerun end-to-end on the user's CDF archive.

## Data acknowledgments

PSP/FIELDS magnetometer data from NASA CDAWeb and the FIELDS team.
Scientific background: [FIELDS data portal](https://fields.ssl.berkeley.edu/data/).
Holdout reference: [Martinović et al. catalog](https://doi.org/10.5281/zenodo.4579821).
