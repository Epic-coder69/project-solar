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
| `magnetic_data_reader.py` | One magnetic-field reader for one-minute and four-samples-per-cycle FIELDS data |
| `magnetic_data_visualization.py` | Daily magnetic-field overview plots |
| `magnetic_data_analysis.py` | Candidate switchback detector; local or exploratory SPC reference |
| `config/survey_selected.csv` | Frozen set of 20 exploratory November 2018 seed windows |
| `results/` | Small set of retained detector catalogs and final validation summary |
| `figures/` | Illustrations and figures used in the research document |

**Experimental scripts, raw comparisons, intermediate results and review images
were removed from `main` to reduce clutter. Nothing was destroyed:** the complete
pre-cleanup snapshot remains on
[`archive/pre-cleanup-2026-10-25`](https://github.com/Epic-coder69/project-solar/tree/archive/pre-cleanup-2026-10-25).

## Setup

Python 3.10+ recommended.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Raw NASA CDF data belong in the git-ignored `data/` cache. The loaders use
NASA CDAWeb and fetch missing files when needed; FIELDS high-resolution data
can be large.

## Basic magnetic-field plots

```bash
python magnetic_data_reader.py --start 2018-11-06 --days 1 --resolution 1min
python magnetic_data_reader.py --start 2018-11-06 --days 1 --resolution 4cyc
python magnetic_data_visualization.py --start 2018-11-06 --days 1
```

Both FIELDS products are read through the same module: `get_day(day, resolution="1min")`
or `get_day(day, resolution="4cyc")`. The detector uses the higher-cadence
measurements for event timing and the one-minute measurements for its slowly
varying magnetic reference. This does **not** manufacture high-resolution
data from one-minute averages.

## Candidate switchback detector

```bash
python magnetic_data_analysis.py --self-test
python magnetic_data_analysis.py --top 1 --reference local --no-plots
```

The full frozen selection is available with `--top 0`; `--top 1` is just a
quick functional run. Results are generated under `data/` and intentionally
are **not automatically committed**.

An **experimental** Parker-spiral comparison is available via
`--reference parker-spc` (requires SWEAP SPC proton velocity/position files).
This is **not** a reproduction of the published SPAN-i-based method.

### Scientific status

- Reading and preprocessing FIELDS magnetic-field vectors: **working**.
- Automatic rotation/switchback *candidate* extraction: **working**.
- Independent validation: **completed exploratory comparison, inconclusive
  detector ranking**; see [holdout findings](results/holdout_findings.md).
- Physically confirmed switchbacks and robust event boundaries: **not yet**.
- Energetic-particle correlation analysis with IS☉IS: **future work**.

**Reproducibility notice:** The high-resolution loading function has been merged into `magnetic_data_reader.py`.
The original local `switchback_survey.py` was not tracked in Git; thus
reproduction of the original event catalogs still needs checking against
actual NASA CDF data on the research Mac. The refactor was tested with
simulated CDF arrays but not yet rerun end-to-end on the user's CDF archive.

## Data acknowledgments

PSP/FIELDS magnetometer data from NASA CDAWeb and the FIELDS team.
Scientific background: [FIELDS data portal](https://fields.ssl.berkeley.edu/data/).
Holdout reference: [Martinović et al. catalog](https://doi.org/10.5281/zenodo.4579821).
