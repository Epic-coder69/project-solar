#!/usr/bin/env python3
"""Project Solar: compare magnetic deflection references on identical PSP data.

Generates overlaid angular diagnostics and an UNLABELED review worksheet.
Never treats detector agreement as confirmation of physical switchbacks.

Requires reference_detector_v2.py, compare_detectors.py, and cached/downloadable
FIELDS 4cyc + one-minute and SWEAP SPC files (using existing project readers).

Published example: Laker et al. 2022, Figure 1, 2019-08-31 ~12:35-12:47 UTC.
https://doi.org/10.1093/mnras/stac2477
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DATA = Path("data")
FIG = Path("figures")
REVIEW = Path("results/reference_review_sheet.csv")


def as_time(value):
    """All project CSV timestamps represent UTC, without timezone offsets."""
    return pd.Timestamp(value).tz_localize(None)


def group_segments(value, date):
    """Decode compare_detectors.py's semicolon-separated UTC clock ranges."""
    if pd.isna(value) or not str(value).strip():
        return []
    intervals = []
    for item in str(value).split(";"):
        if not item.strip():
            continue
        start, end = item.strip().split("–")
        intervals.append((as_time(f"{date}T{start.strip()}"),
                          as_time(f"{date}T{end.strip()}")))
    return intervals


def choose_seed(row, selected):
    """Pick shared seed that generated this group; no time-based invention."""
    date = str(row["date"])
    available = selected[
        selected["start"].astype(str).str.startswith(date)
    ].copy()
    if available.empty:
        raise ValueError(f"No survey seed for {date}")
    midpoint = as_time(row["group_start"]) + (
        as_time(row["group_end"]) - as_time(row["group_start"])
    ) / 2
    delta = (pd.to_datetime(available["start"]) - midpoint).abs()
    return available.loc[delta.idxmin()].to_dict()


def write_review_sheet(rows):
    """Do not overwrite any review decisions already recorded by the user."""
    REVIEW.parent.mkdir(exist_ok=True, parents=True)
    fields = ["group_id", "date", "relationship", "figure",
              "reviewer", "manual_classification", "boundary_quality",
              "notes", "reviewed_utc"]
    fresh = pd.DataFrame(rows)
    for name in fields:
        if name not in fresh:
            fresh[name] = ""
    if REVIEW.exists():
        previous = pd.read_csv(REVIEW, keep_default_na=False)
        known = {str(r["group_id"]): r for r in
                 previous.to_dict("records")}
        for i, row in fresh.iterrows():
            original = known.get(str(row["group_id"]))
            if original:
                for name in ("reviewer", "manual_classification",
                             "boundary_quality", "notes", "reviewed_utc"):
                    fresh.at[i, name] = original.get(name, "")
    fresh[fields].to_csv(REVIEW, index=False)


def diagnostic(row, selected, v2, padding_minutes=3):
    seed = choose_seed(row, selected)
    b, local_ref, mag_ref, _, _ = v2.get_window(seed, reference="local")
    spc_b, spc_ref, _, _, _ = v2.get_window(seed, reference="parker-spc")
    local_angle = v2.deflection(b, local_ref)
    spc_angle = v2.deflection(spc_b, spc_ref).reindex(b.index)

    start = as_time(row["group_start"]) - pd.Timedelta(minutes=padding_minutes)
    end = as_time(row["group_end"]) + pd.Timedelta(minutes=padding_minutes)
    segment = b.loc[start:end]
    la = local_angle.loc[start:end]
    pa = spc_angle.loc[start:end]
    mr = mag_ref.loc[start:end]

    both = la.notna() & pa.notna()
    valid_count = int(both.sum())
    difference = (la[both] - pa[both]).abs()
    median_delta = float(difference.median()) if valid_count else np.nan
    pct_l = float(la[both].ge(45).mean() * 100) if valid_count else np.nan
    pct_p = float(pa[both].ge(45).mean() * 100) if valid_count else np.nan

    fig, axs = plt.subplots(3, 1, figsize=(13, 10),
                            sharex=True, constrained_layout=True)
    for component in ("BR", "BT", "BN"):
        axs[0].plot(segment.index, segment[component], lw=.8,
                    label=component)
    axs[0].legend(loc="upper right", ncol=3)
    axs[0].set_ylabel("RTN field (nT)")

    axs[1].plot(segment.index, segment.Bmag, lw=1,
                color="black", label="Measured |B|")
    axs[1].plot(mr.index, mr, lw=.9, ls="--",
                color="purple", label="Rolling median of scalar |B|")
    axs[1].legend(loc="upper right")
    axs[1].set_ylabel("|B| (nT)")

    axs[2].plot(la.index, la, lw=1.2, label="Local baseline", color="#1674BD")
    axs[2].plot(pa.index, pa, lw=1.1, label="Parker-SPC approx.", color="#CC7A1B")
    axs[2].axhline(30, ls="--", color="gray", lw=.8,
                   label="30° boundary")
    axs[2].axhline(45, ls="--", color="red", lw=.8,
                   label="45° trigger")
    axs[2].set_ylabel("Deflection (degrees)")
    axs[2].set_xlabel("UTC")
    axs[2].set_ylim(-2, 183)

    for i, (field, color, name) in enumerate(
        (("local_intervals", "#1674BD", "Local detections"),
         ("parker_spc_intervals", "#CC7A1B", "SPC detections"))
    ):
        for j, (begin, finish) in enumerate(group_segments(row[field], row["date"])):
            axs[2].axvspan(begin, finish, color=color, alpha=.13,
                           label=name if j == 0 else None)
            # Bottom marker allows viewing boundaries even when spans overlap.
            y = 176 - i * 6
            axs[2].hlines(y, begin, finish, colors=color, lw=4)

    for ax in axs:
        ax.grid(alpha=.2)
        ax.set_xlim(start, end)
    axs[2].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
    axs[2].legend(ncol=3, loc="upper center", fontsize=8)
    fig.suptitle(
        f"Reference-method boundary audit | group {row['group_id']} | "
        f"{row['relationship']} | {row['date']} UTC\n"
        f"Common valid 1-s samples: {valid_count} | "
        f"Median angle difference: {median_delta:.1f}°"
    )

    FIG.mkdir(exist_ok=True, parents=True)
    output = FIG / f"reference_audit_group_{int(row['group_id']):02d}.png"
    fig.savefig(output, dpi=150)
    plt.close(fig)
    return dict(
        group_id=int(row["group_id"]), date=row["date"],
        relationship=row["relationship"], figure=output.as_posix(),
        paired_valid_seconds=valid_count,
        median_absolute_angle_difference_deg=round(median_delta, 2),
        local_fraction_above_45_pct=round(pct_l, 2),
        parker_spc_fraction_above_45_pct=round(pct_p, 2),
    )


def published_example(v2):
    """Independent published illustrative interval; NOT labeled ground truth."""
    seed = dict(start="2019-08-31T12:41:00",
                end="2019-08-31T12:47:00", category="literature_example")
    b, lref, magref, _, _ = v2.get_window(seed, reference="local")
    pb, pref, _, _, _ = v2.get_window(seed, reference="parker-spc")
    a = v2.deflection(b, lref)
    c = v2.deflection(pb, pref).reindex(b.index)
    start, end = pd.Timestamp("2019-08-31T12:35:00"), pd.Timestamp("2019-08-31T12:50:00")
    fig, axes = plt.subplots(3, 1, figsize=(12, 8),
                             sharex=True, constrained_layout=True)
    for comp in ("BR", "BT", "BN"):
        axes[0].plot(b.loc[start:end].index, b.loc[start:end, comp],
                     label=comp, lw=.8)
    axes[0].legend(ncol=3); axes[0].set_ylabel("RTN B (nT)")
    axes[1].plot(b.loc[start:end].index, b.Bmag.loc[start:end],
                 color="black", label="Measured |B|")
    axes[1].plot(magref.loc[start:end].index, magref.loc[start:end],
                 ls="--", color="purple", label="Reference |B|")
    axes[1].legend(); axes[1].set_ylabel("|B| (nT)")
    axes[2].plot(a.loc[start:end].index, a.loc[start:end], label="Local", color="#1674BD")
    axes[2].plot(c.loc[start:end].index, c.loc[start:end], label="Parker-SPC approx.", color="#CC7A1B")
    axes[2].axhline(30, ls="--", color="gray")
    axes[2].axhline(45, ls="--", color="red")
    axes[2].legend(); axes[2].set_ylabel("Angle (degrees)")
    axes[2].set_xlabel("UTC")
    for ax in axes:
        ax.set_xlim(start, end); ax.grid(alpha=.2)
    fig.suptitle("External reference check — Laker et al. Fig. 1, 2019-08-31\n"
                 "Our two approximations (NOT a reproduction of published curve)")
    FIG.mkdir(exist_ok=True, parents=True)
    destination = FIG / "reference_audit_laker2019_example.png"
    fig.savefig(destination, dpi=150); plt.close(fig)
    print("Published-example plot:", destination)


def self_test():
    example = dict(date="2018-11-05",
                   local_intervals="13:40:00–13:42:00; 13:42:20–13:43:00")
    result = group_segments(example["local_intervals"], example["date"])
    assert len(result) == 2
    assert (result[0][1] - result[0][0]).total_seconds() == 120
    assert group_segments("", example["date"]) == []
    sample = pd.DataFrame([
        {"start":"2018-11-05T13:50:30","end":"2018-11-05T13:51:00"},
        {"start":"2018-11-07T03:32:30","end":"2018-11-07T03:35:00"}
    ])
    choice = choose_seed(
        dict(date="2018-11-05", group_start="2018-11-05T13:45:00",
             group_end="2018-11-05T14:00:00"), sample)
    assert choice["start"] == "2018-11-05T13:50:30"
    print("PASS: UTC interval parsing and seed selection")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--all", action="store_true",
                        help="Plot concordant as well as disagreement groups")
    parser.add_argument("--published-example", action="store_true",
                        help="Also inspect Laker et al. Aug 31 2019 example")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return

    group_path = DATA / "reference_overlap_groups.csv"
    seed_path = DATA / "survey_selected.csv"
    if not group_path.exists() or not seed_path.exists():
        parser.error("Run compare_detectors.py first; need overlap groups and survey_selected.csv")
    groups = pd.read_csv(group_path)
    if not args.all:
        groups = groups[groups.relationship != "one_to_one"]
    seeds = pd.read_csv(seed_path)
    import reference_detector_v2 as v2

    outputs = []
    for row in groups.to_dict("records"):
        try:
            output = diagnostic(row, seeds, v2)
            outputs.append(output)
            print(f"[GROUP {row['group_id']}] {row['relationship']}: {output['figure']}")
        except Exception as exc:
            print(f"[FAILED group {row['group_id']}] {exc}")
    if outputs:
        DATA.mkdir(exist_ok=True)
        pd.DataFrame(outputs).to_csv(DATA / "reference_angle_diagnostics.csv", index=False)
        write_review_sheet(outputs)
        print("Wrote data/reference_angle_diagnostics.csv and", REVIEW)
    if args.published_example:
        try:
            published_example(v2)
        except Exception as exc:
            print(f"[Published example unavailable] {exc}")
    if not outputs and not args.published_example:
        raise SystemExit("No diagnostic plots completed.")


if __name__ == "__main__":
    main()
