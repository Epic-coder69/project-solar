"""Generate close-up plots of preliminary candidates."""

import argparse
from pathlib import Path

import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from psp_fields import get_day

parser = argparse.ArgumentParser()
parser.add_argument("--top", type=int, default=5)
args = parser.parse_args()

catalog = Path("data/candidates_1min.csv")

if not catalog.exists():
    raise SystemExit("Run screen_rotations.py first")

candidates = pd.read_csv(catalog)

if candidates.empty:
    raise SystemExit(
        "No candidates found. Try a longer interval."
    )

Path("figures").mkdir(exist_ok=True)

selected = candidates.sort_values(
    "max_angle_deg", ascending=False
).head(args.top)

for _, event in selected.iterrows():
    start = pd.Timestamp(event.start_utc)
    end = pd.Timestamp(event.end_utc)

    left = start - pd.Timedelta(minutes=20)
    right = end + pd.Timedelta(minutes=20)

    day = start.date()
    df = get_day(day)

    if left.date() < day:
        df = pd.concat([
            get_day(day - pd.Timedelta(days=1)), df
        ])

    if right.date() > day:
        df = pd.concat([
            df, get_day(day + pd.Timedelta(days=1))
        ])

    segment = df.loc[left:right]

    fig, axes = plt.subplots(
        2, 1, figsize=(12, 6),
        sharex=True, constrained_layout=True
    )

    for component in ["BR", "BT", "BN"]:
        axes[0].plot(
            segment.index,
            segment[component],
            label=component,
            lw=1
        )

    axes[0].legend(ncol=3)
    axes[0].set_ylabel("Magnetic components (nT)")

    axes[1].plot(
        segment.index,
        segment.Bmag,
        color="black"
    )
    axes[1].set_ylabel("|B| (nT)")
    axes[1].set_xlabel("UTC")

    for ax in axes:
        ax.grid(alpha=0.3)
        ax.axvspan(
            start, end,
            color="orange", alpha=0.2
        )

    fig.suptitle(
        f"Magnetic-rotation candidate: "
        f"{start:%Y-%m-%d %H:%M} UTC"
    )

    output = (
        Path("figures") /
        f"candidate_{start:%Y%m%d_%H%M}.png"
    )

    fig.savefig(output, dpi=160)
    plt.close(fig)
    print("Saved:", output)
