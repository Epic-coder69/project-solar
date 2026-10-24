#!/usr/bin/env python3
"""Independent switchback-catalog comparison (not a ground-truth accuracy test).

Martinovic et al. (2021), Zenodo DOI 10.5281/zenodo.4579821.
Compare independently selected switchbacks against local and Parker-SPC
candidate detectors over identical, successfully evaluated observation windows.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

import h5py
import numpy as np
import pandas as pd

ZENODO = "https://zenodo.org/records/4579821/files/ApJ_SB_Database.h5?download=1"
MD5 = "5046a6964a2f17383c243892202ad5c5"


def download(path):
    import requests
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and hashlib.md5(path.read_bytes()).hexdigest() == MD5:
        return
    response = requests.get(ZENODO, timeout=120)
    response.raise_for_status()
    path.write_bytes(response.content)
    digest = hashlib.md5(path.read_bytes()).hexdigest()
    if digest != MD5:
        path.unlink(missing_ok=True)
        raise ValueError(f"Zenodo checksum mismatch: {digest}; expected {MD5}")


def inspect_h5(path):
    with h5py.File(path) as f:
        def show(name, item):
            if isinstance(item, h5py.Dataset):
                preview = np.asarray(item[0]).tolist() if item.ndim and len(item) else None
                print(f"DATASET {name}: shape={item.shape}, dtype={item.dtype}, first={str(preview)[:230]}")
            else:
                print("GROUP", name)
        f.visititems(show)


def read_table(path, key):
    with h5py.File(path) as f:
        if key not in f:
            raise ValueError(f"HDF5 path not found: {key}. Run --inspect.")
        item = f[key]
        if isinstance(item, h5py.Dataset):
            array = item[()]
            if array.dtype.names:
                return pd.DataFrame.from_records(array)
            if array.ndim == 2:
                return pd.DataFrame(array, columns=[str(i) for i in range(array.shape[1])])
            raise ValueError(f"Unsupported shape {array.shape} at {key}; select a table")
        cols = {name: x[()] for name, x in item.items()
                if isinstance(x, h5py.Dataset) and x.ndim == 1}
        if not cols or len({len(v) for v in cols.values()}) != 1:
            raise ValueError(f"{key} is not a group of equal-length columns; run --inspect")
        return pd.DataFrame(cols)


def pick_fields(frame, start, end):
    fields = {re.sub("[^a-z0-9]", "", str(k).lower()): str(k) for k in frame.columns}
    def resolve(given, options, side):
        if given is not None:
            if given not in frame.columns:
                raise ValueError(f"{side} field {given!r} absent from {list(frame.columns)}")
            return given
        matches = set(fields[k] for k in options if k in fields)
        if len(matches) != 1:
            raise ValueError(f"Unclear {side} field among {list(frame.columns)}. "
                             f"Supply --{side}-field explicitly.")
        return matches.pop()
    return (
        resolve(start, ("sbstart", "switchbackstart", "startutc", "starttime", "tstart"), "start"),
        resolve(end, ("sbend", "switchbackend", "endutc", "endtime", "tend"), "end"),
    )


def utc_times(values, unit):
    s = pd.Series(values).map(lambda x: x.decode("utf-8") if isinstance(x, bytes) else x)
    if unit == "iso":
        return pd.to_datetime(s, utc=True, errors="coerce").dt.tz_localize(None)
    if unit in ("unix_s", "unix_ms", "unix_ns"):
        u = {"unix_s": "s", "unix_ms": "ms", "unix_ns": "ns"}[unit]
        return pd.to_datetime(s, utc=True, unit=u, errors="coerce").dt.tz_localize(None)
    if unit == "mjd":
        return pd.Timestamp("1858-11-17") + pd.to_timedelta(pd.to_numeric(s), unit="D")
    if unit == "jd":
        return pd.Timestamp("1858-11-17") + pd.to_timedelta(pd.to_numeric(s)-2400000.5, unit="D")
    raise ValueError(f"Unrecognized time unit {unit}")


def source_events(args):
    if args.catalog_csv:
        frame = pd.read_csv(args.catalog_csv)
    else:
        path = Path(args.h5)
        if not args.no_download:
            download(path)
        if args.inspect:
            inspect_h5(path)
            return None
        if not args.h5_key:
            inspect_h5(path)
            raise ValueError("Specify --h5-key and event time fields from inspected HDF5; do not guess.")
        frame = read_table(path, args.h5_key)
    if args.inspect:
        print("CSV columns:", list(frame.columns))
        print(frame.head(4).to_string())
        return None
    first, last = pick_fields(frame, args.start_field, args.end_field)
    begin, finish = utc_times(frame[first], args.time_unit), utc_times(frame[last], args.time_unit)
    result = pd.DataFrame({
        "reference_id": [f"M{i+1:04d}" for i in range(len(frame))],
        "start": begin, "end": finish,
    })
    invalid = result.start.isna() | result.end.isna() | (result.end <= result.start)
    if invalid.any():
        raise ValueError(f"{invalid.sum()} invalid timestamps; check time fields and --time-unit")
    return result.sort_values("start").reset_index(drop=True)


def merge_intervals(windows):
    merged = []
    for a, b in sorted(windows):
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    return merged


def contained(start, end, windows):
    return any(a <= start and end <= b for a, b in windows)


def evaluated_windows(seeds, audited_seeds, minutes):
    matched = seeds[seeds.start.astype(str).isin(audited_seeds)].copy()
    if len(matched) != len(audited_seeds):
        raise ValueError("Seed start/end metadata missing; need original data/survey_selected.csv")
    begin = pd.to_datetime(matched.start, utc=True).dt.tz_localize(None)
    finish = pd.to_datetime(matched.end.fillna(matched.start), utc=True).dt.tz_localize(None)
    finish = pd.concat([begin, finish], axis=1).max(axis=1)
    pad = pd.Timedelta(minutes=minutes)
    return merge_intervals(list(zip(begin-pad, finish+pad)))


def detector_candidates(path, allowed_seeds, windows):
    df = pd.read_csv(path)
    if df.empty:
        return pd.DataFrame(columns=["start", "end", "catalog_id"])
    if "source_seed_utcs" not in df:
        df["source_seed_utcs"] = df.seed_utc
    df["source_seed_utcs"] = df.source_seed_utcs.fillna(df.seed_utc)
    df = df[df.source_seed_utcs.apply(
        lambda value: bool(set(str(value).split(";")) & allowed_seeds))].copy()
    df["start"] = pd.to_datetime(df.event_start_utc, utc=True).dt.tz_localize(None)
    df["end"] = pd.to_datetime(df.event_end_utc, utc=True).dt.tz_localize(None)
    df = df[df.apply(lambda r: contained(r.start, r.end, windows), axis=1)].copy()
    df["catalog_id"] = [f"D{i+1:04d}" for i in range(len(df))]
    return df.reset_index(drop=True)


def overlap(a, b, c, d):
    return max(0., (min(b, d)-max(a, c)).total_seconds())


def match_events(refs, detections, threshold):
    """Maximum-cardinality one-to-one event assignment; no claim of ground truth."""
    edges = {}
    for i, r in refs.iterrows():
        options = []
        for j, d in detections.iterrows():
            seconds = overlap(r.start, r.end, d.start, d.end)
            shorter = min((r.end-r.start).total_seconds(), (d.end-d.start).total_seconds())
            if shorter > 0 and seconds/shorter >= threshold:
                options.append((j, seconds))
        edges[i] = [j for j, _ in sorted(options, key=lambda x: -x[1])]
    owner = {}
    def augment(i, used):
        for j in edges[i]:
            if j in used:
                continue
            used.add(j)
            if j not in owner or augment(owner[j], used):
                owner[j] = i
                return True
        return False
    for i in sorted(edges, key=lambda x: len(edges[x])):
        augment(i, set())
    return [(i, j, overlap(refs.loc[i, "start"], refs.loc[i, "end"],
                           detections.loc[j, "start"], detections.loc[j, "end"]))
            for j, i in sorted(owner.items())]


def self_test():
    t = pd.Timestamp("2018-11-01")
    td = pd.Timedelta
    r = pd.DataFrame([{"start": t, "end": t+td(seconds=100)},
                      {"start": t+td(seconds=200), "end": t+td(seconds=300)}])
    d = pd.DataFrame([{"start": t-td(seconds=10), "end": t+td(seconds=90)},
                      {"start": t+td(seconds=190), "end": t+td(seconds=290)}])
    assert len(match_events(r, d, .5)) == 2
    assert len(match_events(r, d, .95)) == 0
    assert merge_intervals([(t, t+td(seconds=30)),
                            (t+td(seconds=20), t+td(seconds=40))]) == [(t, t+td(seconds=40))]
    assert not contained(t, t+td(seconds=50), [(t, t+td(seconds=10))])
    assert pick_fields(pd.DataFrame(columns=["SB_Start", "SB_End"]), None, None) == ("SB_Start", "SB_End")
    print("PASS: overlap thresholds, one-to-one matching, exposure merging, time-field detection")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--h5", default="data/external/ApJ_SB_Database.h5")
    p.add_argument("--catalog-csv", help="Normalized or source CSV with event boundary fields")
    p.add_argument("--h5-key", help="HDF5 dataset/group path determined using --inspect")
    p.add_argument("--start-field")
    p.add_argument("--end-field")
    p.add_argument("--time-unit", choices=("iso", "unix_s", "unix_ms", "unix_ns", "mjd", "jd"),
                   default="iso")
    p.add_argument("--inspect", action="store_true")
    p.add_argument("--no-download", action="store_true")
    p.add_argument("--scope", choices=("shared", "local"), default="shared")
    p.add_argument("--threshold", type=float, default=.5)
    p.add_argument("--padding-minutes", type=float, default=40)
    p.add_argument("--seed-file", default="data/survey_selected.csv")
    p.add_argument("--results-dir", default="results")
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    if args.self_test:
        return self_test()
    if not (0 < args.threshold <= 1):
        p.error("--threshold must be in (0,1]")
    reference = source_events(args)
    if reference is None:
        return
    out = Path(args.results_dir)
    out.mkdir(parents=True, exist_ok=True)
    seeds_path = Path(args.seed_file)
    if not seeds_path.exists():
        p.error(f"Missing {seeds_path}. Rerun the 30-day survey.")
    seeds = pd.read_csv(seeds_path)
    methods = ("local", "parker_spc") if args.scope == "shared" else ("local",)
    events = {k: out/f"corrected_{k}_events.csv" for k in methods}
    audits = {k: out/f"corrected_{k}_audit.csv" for k in methods}
    for path in list(events.values()) + list(audits.values()):
        if not path.exists():
            p.error(f"Missing {path}; copy your detector outputs to results/")
    audited = {}
    for k in methods:
        df = pd.read_csv(audits[k])
        audited[k] = set(df.loc[df.error.fillna("").eq(""), "seed_utc"].astype(str))
    shared = set.intersection(*(audited[k] for k in methods))
    if not shared:
        p.error("No common successfully evaluated seeds")
    windows = evaluated_windows(seeds, shared, args.padding_minutes)
    eligible = reference[reference.apply(
        lambda r: contained(r.start, r.end, windows), axis=1)].copy().reset_index(drop=True)
    eligible.to_csv(out/"external_reference_eligible.csv", index=False)
    rows, summary = [], []
    for k in methods:
        det = detector_candidates(events[k], shared, windows)
        matched = match_events(eligible, det, args.threshold)
        for i, j, seconds in matched:
            r, d = eligible.loc[i], det.loc[j]
            rows.append(dict(
                method=k, reference_id=r.reference_id,
                reference_start=r.start.isoformat(), reference_end=r.end.isoformat(),
                detection_start=d.start.isoformat(), detection_end=d.end.isoformat(),
                overlap_seconds=round(seconds, 1),
                reference_overlap_fraction=round(seconds/(r.end-r.start).total_seconds(), 3),
                detection_overlap_fraction=round(seconds/(d.end-d.start).total_seconds(), 3)))
        summary.append(dict(
            method=k, shared_seed_windows=len(shared), merged_exposure_windows=len(windows),
            external_catalog_total=len(reference), external_events_fully_in_scope=len(eligible),
            detector_intervals_in_scope=len(det), matched_external_events=len(matched),
            unmatched_external_events=len(eligible)-len(matched),
            unmatched_detector_intervals=len(det)-len(matched),
            external_catalog_recovery_fraction=round(len(matched)/len(eligible), 4)
                if len(eligible) else None,
            candidate_catalog_agreement_fraction=round(len(matched)/len(det), 4)
                if len(det) else None))
    pd.DataFrame(rows, columns=[
        "method", "reference_id", "reference_start", "reference_end",
        "detection_start", "detection_end", "overlap_seconds",
        "reference_overlap_fraction", "detection_overlap_fraction"
    ]).to_csv(out/"external_catalog_matches.csv", index=False)
    pd.DataFrame(summary).to_csv(out/"external_catalog_summary.csv", index=False)
    (out/"external_catalog_summary.json").write_text(json.dumps({
        "source": "Martinovic et al. 2021 DOI:10.5281/zenodo.4579821",
        "scope": args.scope, "threshold": args.threshold,
        "criterion": "fraction of shorter temporal interval",
        "exposure_windows": [(a.isoformat(), b.isoformat()) for a, b in windows],
        "results": summary}, indent=2)+"\n")
    print(pd.DataFrame(summary).to_string(index=False))
    print("Saved: external_catalog_matches.csv, external_catalog_summary.csv/.json, external_reference_eligible.csv")
    print("WARNING: Positive-only visually selected catalog; agreement is NOT detector sensitivity, precision, specificity, or false-positive rate.")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError) as error:
        sys.exit(f"ERROR: {error}")
