#!/usr/bin/env python3
"""Reassess catalog agreement under several temporal-matching criteria.

Inputs are the already-committed output CSVs from external_switchback_validation.py.
Uses EXACTLY the recorded shared observation windows, with one-to-one maximum-
cardinality matching and no assumption that catalog positives are exhaustive.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


METHODS = ("local", "parker_spc")
RULES = (
    ("shorter_overlap_50pct", "Overlap >= 50% of shorter interval"),
    ("both_overlap_25pct", "Reference >= 50% AND detection >= 25%"),
    ("iou_10pct", "Intersection/union >= 10%"),
    ("iou_25pct", "Intersection/union >= 25%"),
    ("iou_50pct", "Intersection/union >= 50%"),
    ("both_overlap_30sec", "Overlap >= 30 s; reference >= 50%, detection >= 25%"),
)


def timestamp(x):
    return pd.Timestamp(x, tz="UTC") if pd.Timestamp(x).tz is None else pd.Timestamp(x).tz_convert("UTC")


def in_exposure(start, end, windows):
    return any(a <= start and end <= b for a, b in windows)


def load_data(directory):
    directory = Path(directory)
    meta = json.loads((directory / "external_catalog_summary.json").read_text())
    if meta.get("scope") != "shared":
        raise ValueError("Sensitivity comparison requires --scope shared original results; rerun validation.")
    windows = [(timestamp(a), timestamp(b)) for a, b in meta["exposure_windows"]]
    refs = pd.read_csv(directory / "external_reference_eligible.csv")
    refs["start"] = pd.to_datetime(refs.start, utc=True)
    refs["end"] = pd.to_datetime(refs.end, utc=True)
    refs = refs.loc[refs.apply(lambda r: in_exposure(r.start, r.end, windows), axis=1)]
    if not refs.reference_id.is_unique:
        raise ValueError("Duplicate published event IDs in eligible catalog")
    shared_seeds = None
    for name in METHODS:
        audit = pd.read_csv(directory / f"corrected_{name}_audit.csv")
        successful = set(audit.loc[audit.error.fillna("").eq(""), "seed_utc"].astype(str))
        shared_seeds = successful if shared_seeds is None else shared_seeds & successful
    if not shared_seeds:
        raise ValueError("No successfully evaluated shared seed windows")
    records = {}
    for name in METHODS:
        det = pd.read_csv(directory / f"corrected_{name}_events.csv")
        if "source_seed_utcs" not in det:
            det["source_seed_utcs"] = det.seed_utc
        det["source_seed_utcs"] = det.source_seed_utcs.fillna(det.seed_utc)
        det = det.loc[det.source_seed_utcs.apply(
            lambda s: bool(set(str(s).split(";")) & shared_seeds))].copy()
        det["start"] = pd.to_datetime(det.event_start_utc, utc=True)
        det["end"] = pd.to_datetime(det.event_end_utc, utc=True)
        det = det.loc[det.apply(lambda r: in_exposure(r.start, r.end, windows), axis=1)]
        records[name] = det.reset_index(drop=True)
    return refs.reset_index(drop=True), records, windows


def edges_and_matches(refs, det, rule):
    edges = {}
    for i, r in refs.iterrows():
        options = []
        for j, d in det.iterrows():
            overlap = max(0., (min(r.end, d.end) - max(r.start, d.start)).total_seconds())
            rlen = (r.end - r.start).total_seconds()
            dlen = (d.end - d.start).total_seconds()
            if overlap <= 0 or min(rlen, dlen) <= 0:
                continue
            entry = dict(
                ref_index=i, det_index=j,
                reference_id=r.reference_id,
                start_utc=r.start.isoformat(), end_utc=r.end.isoformat(),
                detection_start_utc=d.start.isoformat(), detection_end_utc=d.end.isoformat(),
                overlap_seconds=overlap, shorter_fraction=overlap/min(rlen, dlen),
                reference_fraction=overlap/rlen,
                detection_fraction=overlap/dlen,
                iou=overlap/(rlen+dlen-overlap))
            if qualifies(entry, rule):
                options.append(entry)
        edges[i] = sorted(options, key=lambda e: (-e["iou"], -e["overlap_seconds"],
                                                   str(e["detection_start_utc"])))
    owner = {}

    def assign(i, seen):
        for e in edges[i]:
            j = e["det_index"]
            if j in seen:
                continue
            seen.add(j)
            if j not in owner or assign(owner[j]["ref_index"], seen):
                owner[j] = e
                return True
        return False

    for i in sorted(edges, key=lambda i: (len(edges[i]), i)):
        assign(i, set())
    return sorted(owner.values(), key=lambda e: (str(e["start_utc"]), e["det_index"]))


def qualifies(e, rule):
    if rule == "shorter_overlap_50pct":
        return e["shorter_fraction"] >= .5
    if rule == "both_overlap_25pct":
        return e["reference_fraction"] >= .5 and e["detection_fraction"] >= .25
    if rule == "iou_10pct":
        return e["iou"] >= .1
    if rule == "iou_25pct":
        return e["iou"] >= .25
    if rule == "iou_50pct":
        return e["iou"] >= .5
    if rule == "both_overlap_30sec":
        return (e["overlap_seconds"] >= 30 and e["reference_fraction"] >= .5
                and e["detection_fraction"] >= .25)
    raise ValueError(f"Unknown rule: {rule}")


def self_test():
    t = pd.Timestamp("2018-11-01", tz="UTC")
    ref = pd.DataFrame([
        dict(reference_id="short", start=t, end=t + pd.Timedelta(seconds=4)),
        dict(reference_id="long", start=t + pd.Timedelta(seconds=10),
             end=t + pd.Timedelta(seconds=70))])
    det = pd.DataFrame([
        dict(start=t, end=t + pd.Timedelta(minutes=10)),
        dict(start=t + pd.Timedelta(seconds=10),
             end=t + pd.Timedelta(seconds=70))])
    loose = edges_and_matches(ref, det, "shorter_overlap_50pct")
    strict = edges_and_matches(ref, det, "iou_25pct")
    assert len(loose) == 2
    assert len(strict) == 1 and strict[0]["reference_id"] == "long"
    assert strict[0]["iou"] == 1
    print("PASS: one-to-one matching and long-interval false agreement rejection")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--results-dir", default="results")
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    if args.self_test:
        return self_test()
    output = Path(args.results_dir)
    refs, records, windows = load_data(output)
    summary, details = [], []
    for rule, description in RULES:
        for method in METHODS:
            pairs = edges_and_matches(refs, records[method], rule)
            summary.append(dict(
                rule=rule, description=description, method=method,
                shared_exposure_windows=len(windows),
                published_events=len(refs), candidate_intervals=len(records[method]),
                matched_published_events=len(pairs),
                published_match_fraction=round(len(pairs)/len(refs), 4) if len(refs) else np.nan,
                matched_candidate_fraction=round(len(pairs)/len(records[method]), 4)
                    if len(records[method]) else np.nan))
            details.extend(dict(rule=rule, method=method, **pair) for pair in pairs)
    result = pd.DataFrame(summary)
    result.to_csv(output / "external_match_sensitivity.csv", index=False)
    pd.DataFrame(details).to_csv(output / "external_match_sensitivity_pairs.csv", index=False)
    print(result[["rule", "method", "published_events", "matched_published_events",
                  "published_match_fraction"]].to_string(index=False))
    print("Saved: results/external_match_sensitivity.csv and _pairs.csv")
    print("WARNING: Catalog matching is not an unbiased estimate of recall or precision.")


if __name__ == "__main__":
    main()
