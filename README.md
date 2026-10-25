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
| `psp_fields.py` | Download/read one-minute RTN magnetic fields |
| `psp_highres.py` | Download/read four-samples-per-cycle RTN magnetic fields |
| `overview.py` | Daily magnetic-field overview plots |
| `reference_detector_v2.py` | Candidate switchback detector; local or exploratory SPC reference |
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
python psp_fields.py --start 2018-11-06 --days 1
python overview.py --start 2018-11-06 --days 1
```

## Candidate switchback detector

```bash
python reference_detector_v2.py --self-test
python reference_detector_v2.py --top 1 --reference local --no-plots
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

**Reproducibility notice:** `psp_highres.py` was consolidated during cleanup
because the earlier detector depended on an untracked local
`switchback_survey.py`. The replacement is syntactically checked but must be
tested against the local NASA CDF cache before claiming that rerunning the
pipeline reproduces the previous catalogs exactly.

## Data acknowledgments

PSP/FIELDS magnetometer data from NASA CDAWeb and the FIELDS team.
Scientific background: [FIELDS data portal](https://fields.ssl.berkeley.edu/data/).
Holdout reference: [Martinović et al. catalog](https://doi.org/10.5281/zenodo.4579821).
