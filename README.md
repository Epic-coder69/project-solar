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

### Remove obsolete local folders and logs

Previous Git renames (`M_FIELD` → `FIELDS` and `IS_O_IS` →
`IS⊙IS`) may leave old **local-only** folders containing Python caches or
duplicate files. They no longer exist as separate tracked folders on GitHub.
The root `logs/` directory is also local-only (Git-ignored).

```bash
python project.py cleanup          # Preview every proposed deletion
python project.py cleanup --apply  # Remove only safe duplicates and local logs
```

The tool removes a legacy script **only if byte-for-byte identical** to the
matching current script; distinct or modified files are preserved and
reported. It also discards rebuildable caches inside legacy folders and
empty legacy directories, and removes root `logs/` and root `*.log` files.
**It never touches `data/`**, the 599-event catalog, NASA CDF/NPZ files,
canonical source code or scientific results. If a legacy folder remains, it
contains files the tool intentionally refused to delete. Check its warnings.

### Most useful commands

```bash
# Run every offline instrument self-test (no NASA downloads)
python project.py test

# Recalculate coverage from the data already saved locally
python project.py coverage

# View all reader options (does not download anything)
python project.py highres --help
python project.py sweap --help

# Inspect which SWEAP products the reader knows about
python project.py sweap --list-products
```

For *intentional* new EPI-Hi hourly downloads to fill remaining source gaps,
specify products and the event catalog explicitly:

```bash
python project.py highres \
  --events-csv data/full_analysis_nov_dec_2018_candidates_strong_events.csv \
  --products let1-rates3600 het-rates3600 \
  --download --resume
```

Commands after `python project.py` forward arguments to the existing
instrument-specific scripts; they do **not** silently change detector
measurements or integration times.

## Repository guide

```text
project-solar/
├── project.py                    # Stable, ASCII one-command entry point
├── cleanup_local.py              # Safe removal of duplicate folders/logs
├── README.md                     # Current overview (this file)
├── The_Research.md               # Scientific objectives and research notes
├── .gitignore                    # Keeps NASA CDFs, caches and outputs local
├── docs/
│   ├── REPOSITORY_MAP.md         # Every script, ownership, local datasets
│   └── WORKFLOWS.md              # Preserved long-form technical walkthrough
└── code/
    ├── FIELDS/                   # Magnetic-field data and event detection
    ├── IS⊙IS/                    # Particle readers and coverage audits
    └── SWEAP/                    # Exploratory plasma data inspection
```

The instrument folder is literally `code/IS⊙IS/` on GitHub. Older local
commands using `code/IS_O_IS/` refer to the **previous folder name**.
Prefer `project.py` so older path examples do not break again.

- [Repository map and local data policy](docs/REPOSITORY_MAP.md)
- [Detailed acquisition and analysis workflows](docs/WORKFLOWS.md)
- [Research objectives and scientific background](The_Research.md)

## Data preservation and reproducibility

**`data/` is ignored by Git.** Your 599-event CSV, original NASA CDF files,
compressed per-day particle partitions and coverage report CSV/JSONs remain
on the analysis computer. A `git pull` does not fetch or delete them.
Keep an independent backup of this directory before any restructuring.
Never run `git clean -fdx` or force-reset your working tree to tidy folders.

Small synthetic tests run automatically in **GitHub Actions** on code
changes. Those checks validate software behavior, **not** scientific
calibration, true data availability or causal particle acceleration.

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
