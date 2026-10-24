"""Preliminary magnetic-rotation screening."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from psp_fields import days, get_day

parser = argparse.ArgumentParser()
parser.add_argument("--start", default="2018-11-06")
parser.add_argument("--days", type=int, default=1)
args = parser.parse_args()

rows = []

for day in days(args.start, args.days):
    df = get_day(day)
    B = df[["BR", "BT", "BN"]]

    # Local magnetic-field reference: 61-minute median
    baseline = B.rolling(
        61, center=True, min_periods=30
    ).median()

    baseline_mag = np.sqrt(
        (baseline ** 2).sum(axis=1, min_count=3)
    )

    cosine = (
        (B * baseline).sum(axis=1, min_count=3)
        / (df.Bmag * baseline_mag)
    )

    angle = np.degrees(
        np.arccos(cosine.clip(-1, 1))
    )

    typical_mag = df.Bmag.rolling(
        61, center=True, min_periods=30
    ).median()

    relative_change = (
        (df.Bmag - typical_mag).abs()
        / typical_mag
    )

    # Preliminary criteria:
    # rotation > 100 degrees
    # radial component reverses relative to local reference
    # reference field has meaningful radial component
    # field magnitude differs < 25% from local median
    mask = (
        (angle > 100)
        & (df.BR * baseline.BR < 0)
        & (baseline.BR.abs() > 0.25 * baseline_mag)
        & (relative_change < 0.25)
    )

    positions = np.flatnonzero(
        mask.fillna(False).to_numpy()
    )

    if not len(positions):
        print(day, "No candidates passing screening")
        continue

    groups = np.split(
        positions,
        np.where(np.diff(positions) > 2)[0] + 1
    )

    for group in groups:
        if len(group) < 2:
            continue

        i, j = group[0], group[-1]

        rows.append({
            "start_utc": df.index[i].isoformat(),
            "end_utc": df.index[j].isoformat(),
            "flagged_minutes": len(group),
            "max_angle_deg": round(
                float(angle.iloc[group].max()), 1
            ),
            "max_relative_magnitude_change": round(
                float(relative_change.iloc[group].max()), 3
            )
        })

    print(day, "Groups examined:", len(groups))

output = Path("data/candidates_1min.csv")
output.parent.mkdir(exist_ok=True)

columns = [
    "start_utc",
    "end_utc",
    "flagged_minutes",
    "max_angle_deg",
    "max_relative_magnitude_change"
]

pd.DataFrame(
    rows, columns=columns
).to_csv(output, index=False)

print("Saved", len(rows), "PRELIMINARY candidates")
print("Catalog:", output)
