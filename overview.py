"""Plot daily FIELDS magnetic-field time series."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from psp_fields import days, get_day

parser = argparse.ArgumentParser()
parser.add_argument("--start", default="2018-11-06")
parser.add_argument("--days", type=int, default=1)
args = parser.parse_args()

Path("figures").mkdir(exist_ok=True)

for day in days(args.start, args.days):
    df = get_day(day)

    fig, axes = plt.subplots(
        4, 1, figsize=(13, 9),
        sharex=True, constrained_layout=True
    )

    for ax, name in zip(
        axes, ["BR", "BT", "BN", "Bmag"]
    ):
        ax.plot(df.index, df[name], lw=0.7)
        ax.set_ylabel(f"{name} (nT)")
        ax.grid(alpha=0.3)

        if name != "Bmag":
            ax.axhline(
                0, color="gray", ls="--", lw=0.7
            )

    axes[0].set_title(
        f"PSP FIELDS MAG (RTN) — {day} — 1-minute data"
    )
    axes[-1].set_xlabel("UTC")

    output = (
        Path("figures") /
        f"fields_{day:%Y%m%d}_1min.png"
    )

    fig.savefig(output, dpi=160)
    plt.close(fig)
    print("Saved:", output)
