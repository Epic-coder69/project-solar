"""
PROJECT SOLAR — Literature-threshold switchback benchmark

Based on event-selection thresholds used by Laker et al. (2022):
45-degree trigger, 30-degree boundaries,
10-second minimum duration, merge gaps below 20 seconds.

IMPORTANT:
This implementation uses a 3-hour local magnetic-field median
instead of the SWEAP/ephemeris-derived Parker spiral.
It is a benchmark, NOT an exact replication of Laker et al.
"""

import argparse
from datetime import timedelta
from functools import lru_cache
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from psp_fields import get_day

try:
    from switchback_survey import read_high
except ImportError as exc:
    raise SystemExit(
        "Missing switchback_survey.py in the project folder."
    ) from exc


DATA = Path("data")
FIG = Path("figures")

TRIGGER = 45.0
EDGE = 30.0
MIN_SECONDS = 10
MERGE_SECONDS = 20


@lru_cache(maxsize=4)
def high_for_day(day):
    return read_high(day)


@lru_cache(maxsize=5)
def background_for_day(day):
    blocks = []

    for d in (
        day - timedelta(days=1),
        day,
        day + timedelta(days=1)
    ):
        try:
            blocks.append(
                get_day(d)[["BR", "BT", "BN"]]
            )
        except Exception as exc:
            print(f"[REFERENCE] Missing {d}: {exc}")

    if not blocks:
        raise RuntimeError(
            "No one-minute data for reference direction"
        )

    b = pd.concat(blocks).sort_index()
    b = b[~b.index.duplicated(keep="last")]

    return b.rolling(
        "181min",
        center=True,
        min_periods=30
    ).median()


def angle_series(df, bg):
    field = df[["BR", "BT", "BN"]].to_numpy(float)
    ref = bg[["BR", "BT", "BN"]].to_numpy(float)

    norms = (
        np.linalg.norm(field, axis=1)
        * np.linalg.norm(ref, axis=1)
    )

    cosine = np.full(len(df), np.nan)

    valid = (
        np.isfinite(norms)
        & (norms > 1e-8)
        & np.isfinite(field).all(axis=1)
        & np.isfinite(ref).all(axis=1)
    )

    cosine[valid] = (
        np.sum(field[valid] * ref[valid], axis=1)
        / norms[valid]
    )

    return pd.Series(
        np.degrees(
            np.arccos(np.clip(cosine, -1, 1))
        ),
        index=df.index
    )


def find_events(angle):
    hit = (angle >= EDGE).fillna(False).to_numpy(bool)
    good = angle.notna().to_numpy(bool)

    indices = np.flatnonzero(hit)

    if len(indices) == 0:
        return []

    spans = []

    for group in np.split(
        indices,
        np.where(np.diff(indices) > 1)[0] + 1
    ):
        spans.append(
            (int(group[0]), int(group[-1]))
        )

    merged = [spans[0]]

    for begin, finish in spans[1:]:
        old_begin, old_finish = merged[-1]

        gap = (
            angle.index[begin]
            - angle.index[old_finish]
        ).total_seconds()

        if (
            gap < MERGE_SECONDS
            and good[old_finish:begin + 1].all()
        ):
            merged[-1] = (old_begin, finish)
        else:
            merged.append((begin, finish))

    events = []

    for begin, finish in merged:
        duration = (
            angle.index[finish]
            - angle.index[begin]
        ).total_seconds() + 1

        peak = float(
            angle.iloc[begin:finish + 1].max()
        )

        if duration >= MIN_SECONDS and peak >= TRIGGER:
            events.append(
                (begin, finish, duration, peak)
            )

    return events


def make_window(seed, minutes=40):
    center = pd.Timestamp(seed["start"])

    stop = (
        pd.Timestamp(seed["end"])
        if pd.notna(seed.get("end"))
        else center
    )

    if stop < center:
        stop = center

    left = center - pd.Timedelta(minutes=minutes)
    right = stop + pd.Timedelta(minutes=minutes)

    blocks = []

    dates = pd.date_range(
        left.normalize(),
        right.normalize(),
        freq="D"
    )

    for d in dates:
        try:
            blocks.append(high_for_day(d.date()))
        except Exception as exc:
            print(
                f"[HIGH-RES] Missing {d.date()}: {exc}"
            )

    if not blocks:
        raise RuntimeError(
            "No high-resolution data for interval"
        )

    combined = pd.concat(blocks).sort_index()
    combined = combined[
        ~combined.index.duplicated(keep="last")
    ]

    raw = combined.loc[
        left:right,
        ["BR", "BT", "BN"]
    ]

    if len(raw) < 20:
        raise RuntimeError(
            "Insufficient high-resolution samples"
        )

    # One-second median bins; do not interpolate missing data.
    b = raw.resample("1s").median()

    b["Bmag"] = np.sqrt(
        (b[["BR", "BT", "BN"]] ** 2).sum(
            axis=1,
            min_count=3
        )
    )

    ref_blocks = [
        background_for_day(d.date())
        for d in dates
    ]

    ref = pd.concat(ref_blocks).sort_index()
    ref = ref[
        ~ref.index.duplicated(keep="last")
    ]

    ref = ref.reindex(
        b.index,
        method="nearest",
        tolerance=pd.Timedelta(seconds=35)
    )

    angle = angle_series(b, ref)

    return b, ref, angle


def check_seed(seed, number):
    b, bg, angle = make_window(seed)

    groups = find_events(angle)

    center = pd.Timestamp(seed["start"])

    stop = (
        pd.Timestamp(seed["end"])
        if pd.notna(seed.get("end"))
        else center
    )

    category = str(
        seed.get("category", "unknown")
    )

    results = []

    for i, j, duration, peak in groups:
        start = angle.index[i]
        end = angle.index[j]

        # Only associate nearby events with the original seed.
        if (
            end < center - pd.Timedelta(minutes=15)
            or start > stop + pd.Timedelta(minutes=15)
        ):
            continue

        segment = b.iloc[i:j + 1]
        ref_segment = bg.iloc[i:j + 1]

        mag_ref = np.sqrt(
            (
                ref_segment[["BR", "BT", "BN"]] ** 2
            ).sum(axis=1, min_count=3)
        )

        deviation = (
            segment.Bmag / mag_ref - 1
        ).abs()

        reversal = (
            (segment.BR * ref_segment.BR) < 0
        ).any()

        results.append({
            "seed_utc": center.isoformat(),
            "seed_category": category,
            "event_start_utc": start.isoformat(),
            "event_end_utc": end.isoformat(),
            "duration_seconds": round(duration, 1),
            "maximum_deflection_deg": round(peak, 1),
            "median_abs_Bmag_deviation": round(
                float(deviation.median()), 3
            ),
            "max_abs_Bmag_deviation": round(
                float(deviation.max()), 3
            ),
            "radial_reversal": bool(reversal),
            "status": "Unvalidated magnetic candidate"
        })

    FIG.mkdir(exist_ok=True)

    filename = (
        f"benchmark_{category}_"
        f"{center:%Y%m%d_%H%M}_{number}.png"
    )

    fig, axes = plt.subplots(
        3, 1,
        figsize=(12, 9),
        sharex=True,
        constrained_layout=True
    )

    for component in ["BR", "BT", "BN"]:
        axes[0].plot(
            b.index,
            b[component],
            lw=0.7,
            label=component
        )

    axes[0].legend(ncol=3)
    axes[0].set_ylabel("Magnetic field (nT)")

    axes[1].plot(
        b.index,
        b.Bmag,
        color="black",
        lw=0.8
    )
    axes[1].set_ylabel("|B| (nT)")

    axes[2].plot(
        angle.index,
        angle,
        lw=0.8
    )

    axes[2].axhline(
        EDGE,
        color="gray",
        linestyle="--",
        label="30° boundary"
    )

    axes[2].axhline(
        TRIGGER,
        color="red",
        linestyle="--",
        label="45° trigger"
    )

    axes[2].legend(ncol=2)
    axes[2].set_ylabel("Deflection angle (°)")
    axes[2].set_xlabel("Time (UTC)")

    for event in results:
        x1 = pd.Timestamp(event["event_start_utc"])
        x2 = pd.Timestamp(event["event_end_utc"])

        for ax in axes:
            ax.axvspan(
                x1, x2,
                color="orange",
                alpha=0.15
            )

    for ax in axes:
        ax.grid(alpha=0.2)

    fig.suptitle(
        f"Literature-threshold benchmark | "
        f"{category} | {center:%Y-%m-%d %H:%M} UTC"
    )

    output = FIG / filename

    fig.savefig(output, dpi=150)
    plt.close(fig)

    return results, str(output)


def self_test():
    ix = pd.date_range(
        "2018-11-06",
        periods=600,
        freq="s"
    )

    a = pd.Series(
        np.zeros(len(ix)),
        index=ix
    )

    a.iloc[120:151] = 140
    events = find_events(a)

    assert len(events) == 1
    assert events[0][2] == 31

    a.iloc[:] = 0
    a.iloc[120:126] = 140

    assert find_events(a) == []

    a.iloc[:] = 0
    a.iloc[120:133] = 120
    a.iloc[140:153] = 120

    events = find_events(a)

    assert len(events) == 1
    assert events[0][2] == 33

    print(
        "PASS: event duration, short-event rejection, "
        "and gap merging"
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--top",
        type=int,
        default=8,
        help="Number of survey seeds; 0 means all"
    )

    parser.add_argument(
        "--self-test",
        action="store_true"
    )

    args = parser.parse_args()

    if args.self_test:
        self_test()
        return

    path = DATA / "survey_selected.csv"

    if not path.exists():
        parser.error(
            "Missing data/survey_selected.csv. "
            "Run switchback_survey.py first."
        )

    selected = pd.read_csv(path)

    if selected.empty:
        parser.error("No selected survey intervals found")

    selected["priority"] = selected.category.map({
        "strong": 0,
        "ambiguous": 1,
        "control": 2
    }).fillna(3)

    selected = selected.sort_values(
        ["priority", "start"]
    )

    if args.top > 0:
        selected = selected.head(args.top)

    events = []
    audit = []

    for i, seed in enumerate(
        selected.to_dict("records"), 1
    ):
        try:
            found, figure = check_seed(seed, i)

            events.extend(found)

            audit.append({
                "seed_utc": seed["start"],
                "seed_category": seed["category"],
                "events_detected_near_seed": len(found),
                "figure": figure,
                "error": ""
            })

            print(
                f"[{i}/{len(selected)}] "
                f"{seed['start']}: "
                f"{len(found)} magnetic candidates"
            )

        except Exception as exc:
            audit.append({
                "seed_utc": seed["start"],
                "seed_category": seed["category"],
                "events_detected_near_seed": 0,
                "figure": "",
                "error": str(exc)
            })

            print(
                f"[ERROR] {seed['start']}: {exc}"
            )

    columns = [
        "seed_utc",
        "seed_category",
        "event_start_utc",
        "event_end_utc",
        "duration_seconds",
        "maximum_deflection_deg",
        "median_abs_Bmag_deviation",
        "max_abs_Bmag_deviation",
        "radial_reversal",
        "status"
    ]

    DATA.mkdir(exist_ok=True)

    pd.DataFrame(
        events, columns=columns
    ).to_csv(
        DATA / "benchmark_events.csv",
        index=False
    )

    pd.DataFrame(
        audit
    ).to_csv(
        DATA / "benchmark_audit.csv",
        index=False
    )

    print(
        f"DONE: {len(events)} candidate intervals; "
        f"{len(audit)} evaluated seeds"
    )

    print(
        "Results: data/benchmark_events.csv, "
        "data/benchmark_audit.csv, "
        "figures/benchmark_*.png"
    )

    if not any(x["figure"] for x in audit):
        raise SystemExit(
            "No figures generated. "
            "Check benchmark_audit.csv."
        )


if __name__ == "__main__":
    main()
