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


## Long scans with event archiving and automatic raw-cache cleanup

**Use the rolling pipeline instead of a 2-month reader-first command when
disk space matters.** It processes one UTC day at a time:

1. Reader searches one day and saves its candidate CSV.
2. Analysis characterizes all candidates and saves events/registry/audit.
3. Analysis archives **every 1-second candidate-window waveform**, including
   the observed RTN components, 181-minute background, magnitude and angle.
4. The pipeline archives **every detected event** at the original native sample
   cadence, including five minutes before/after each event, as gzip CSV.
   No 90-degree or radial-reversal filter is used for archival.
5. Only after all outputs and archives are checked does the optional cleanup
   remove older full-cadence CDF files from **data/fields_full/**. It retains
   the preceding day's files while processing adjacent days and never deletes
   unrelated product caches or files outside the requested observation period.
   Any day with less than 90% valid one-second coverage is kept intact.

First run its safe, local synthetic tests:

\`\`\`bash
python code/magnetic_data_pipeline.py --self-test
python code/magnetic_data_analysis.py --self-test
\`\`\`

Two-month unattended example (November 1, 2018 through January 1, 2019 UTC;
run from the repository root, on a Mac connected to power with lid open):

\`\`\`bash
mkdir -p logs
nohup caffeinate -is ./.venv/bin/python code/magnetic_data_pipeline.py \
  --start 2018-11-01 --duration 2mo --cleanup-raw \
  > logs/rolling_nov_dec_2018.log 2>&1 < /dev/null &
tail -f logs/rolling_nov_dec_2018.log
\`\`\`

To resume after interruption, rerun the same command with \`--resume\`.
To preview the cleanup without actually removing any CDFs, add
\`--dry-run-cleanup\` (requires \`--cleanup-raw\`). Leave out
\`--cleanup-raw\` to keep the original CDFs. Add \`--plots\` only if PNG
figures are needed for every window; PNGs are off by default to save space.

Outputs are saved by day under \`data/rolling_<start>_<end>/\`, including
daily checkpoints, event gzip files, one-second background archives and an
event archive index. At successful completion the folder also contains
\`all_events.csv\`, \`all_registry.csv\`, \`all_audit.csv\` and
\`all_event_archives.csv\`.

**Scientific and data-retention caveats:** Archived event samples are
quality-filtered RTN magnetic vectors, *not* a lossless original NASA CDF.
They do not preserve the full instrument metadata or raw quality-flag history.
Quiet full-cadence intervals not inside event archives are deliberately
discarded from the local cache when cleanup is enabled; original CDFs remain
available for later download from NASA. A five-minute context is retained
around each magnetic event, plus the 1-second local-background traces used
for detection. This pipeline does not download or delete plasma and energetic
particle products. A long 2-month run still needs ample disk capacity and
may stop safely on missing data or archive verification failure; completed
daily checkpoints remain recoverable.

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
