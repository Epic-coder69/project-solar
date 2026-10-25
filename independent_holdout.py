#!/usr/bin/env python3
"""Frozen, independent PSP switchback holdout; see results/holdout_protocol.md.

Prepare without looking at detector output. Then run BOTH unmodified v2
reference methods at the same held-out intervals. Random catalog-clear
intervals are UNLABELED, not certified negatives.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path("results")
PAD = pd.Timedelta(minutes=40)
METHODS = ("local", "parker-spc")
BEGIN, END = pd.Timestamp("2018-11-01"), pd.Timestamp("2018-12-01")


def time(x):
    return pd.Timestamp(x).tz_localize(None)


def intersect(a, b, c, d):
    return a < d and c < b


def collides(a, b, windows):
    return any(intersect(a, b, c, d) for c, d in windows)


def span_from_seed(row):
    a = time(row["start"])
    b = time(row["end"]) if pd.notna(row.get("end")) else a
    return a-PAD, max(a, b)+PAD


def iou(a, b, c, d):
    x = max(0., (min(b, d)-max(a, c)).total_seconds())
    denom = (b-a).total_seconds()+(d-c).total_seconds()-x
    return x/denom if denom > 0 else 0.


def sample(refs, original, coverage, positives, controls, seed):
    rng = np.random.default_rng(seed)
    old = [span_from_seed(r) for r in original.to_dict("records")]
    occupied = list(old)
    e1 = refs[refs.reference_id.str.startswith("E1-")].copy()
    e1["start"] = pd.to_datetime(e1.start).dt.tz_localize(None)
    e1["end"] = pd.to_datetime(e1.end).dt.tz_localize(None)
    if e1.empty or (e1.end <= e1.start).any():
        raise ValueError("Encounter 1 reference event list empty or invalid")
    domain_a = max(BEGIN, e1.start.min())
    domain_b = min(END, e1.end.max())
    allowed = set(coverage.loc[coverage.valid_minutes.ge(1000), "day"].astype(str).str[:10])
    available = e1[(e1.start >= BEGIN) & (e1.end < END)]
    available = available[available.start.dt.strftime("%Y-%m-%d").isin(allowed)]
    result = []
    if positives:
        for i in rng.permutation(len(available)):
            r = available.iloc[int(i)]
            a, b = time(r.start), time(r.end)
            if collides(a-PAD, b+PAD, occupied):
                continue
            result.append(dict(stratum="published_positive", reference_id=r.reference_id,
                               start=a.isoformat(), end=b.isoformat()))
            occupied.append((a-PAD, b+PAD))
            if sum(row["stratum"] == "published_positive" for row in result) == positives:
                break
    # A catalog-clear *30-minute core* is not a verified negative.
    protected = [(time(r.start)-pd.Timedelta(minutes=10),
                  time(r.end)+pd.Timedelta(minutes=10))
                 for r in e1.itertuples(index=False)]
    for _ in range(50000):
        if sum(row["stratum"] == "catalog_clear_unlabeled" for row in result) >= controls:
            break
        center = domain_a + (domain_b-domain_a)*float(rng.random())
        if center.strftime("%Y-%m-%d") not in allowed:
            continue
        a, b = center-pd.Timedelta(minutes=15), center+pd.Timedelta(minutes=15)
        if center-PAD < domain_a or center+PAD > domain_b:
            continue
        if collides(center-PAD, center+PAD, occupied) or collides(a, b, protected):
            continue
        result.append(dict(stratum="catalog_clear_unlabeled", reference_id="",
                           start=a.isoformat(), end=b.isoformat()))
        occupied.append((center-PAD, center+PAD))
    result.sort(key=lambda r: r["start"])
    for i, r in enumerate(result, 1):
        r["id"] = f"H{i:03d}"
        if r["stratum"] == "published_positive":
            r["seed_start"], r["seed_end"] = r["start"], r["end"]
        else:
            mid = time(r["start"])+(time(r["end"])-time(r["start"]))/2
            r["seed_start"] = r["seed_end"] = mid.isoformat()
    if not result:
        raise RuntimeError("No independent intervals are available")
    return pd.DataFrame(result), (domain_a, domain_b)


def prepare(args):
    from external_switchback_validation import read_martinovic_events, download
    manifest = OUT/"holdout_manifest.csv"
    if manifest.exists():
        raise SystemExit("Existing holdout_manifest.csv is frozen; never resample to improve results.")
    h5 = Path("data/external/ApJ_SB_Database.h5")
    download(h5)
    refs = read_martinovic_events(h5)
    seeds = pd.read_csv("data/survey_selected.csv")
    coverage = pd.read_csv("results/survey_coverage.csv")
    sampled, domain = sample(refs, seeds, coverage, args.positives, args.unlabeled, args.seed)
    OUT.mkdir(exist_ok=True)
    sampled.to_csv(manifest, index=False)
    config = {
        "random_seed": args.seed, "positives_requested": args.positives,
        "unlabeled_requested": args.unlabeled, "strata_obtained": sampled.stratum.value_counts().to_dict(),
        "catalog_md5": hashlib.md5(h5.read_bytes()).hexdigest(),
        "domain": [str(x) for x in domain],
        "manifest_sha256": hashlib.sha256(manifest.read_bytes()).hexdigest(),
        "exclusion": "All original survey seeds expanded by 40 min each side; all holdout exposures non-overlapping",
        "control_label": "UNLABELED, free of catalog positives in 30-min core +/- 10-min buffer",
    }
    (OUT/"holdout_settings.json").write_text(json.dumps(config, indent=2)+"\n")
    print("Frozen independent holdout:", config["strata_obtained"])
    print("Events with unavailable FIELDS/SPC will be marked as failures, not negatives.")


def evaluate(row, method, v2):
    seed = dict(start=row["seed_start"], end=row["seed_end"], category="holdout")
    b, reference, _, _, _ = v2.get_window(seed, reference=method)
    angles = v2.deflection(b, reference)
    a, z = time(row["start"]), time(row["end"])
    qc_a = a-pd.Timedelta(minutes=5) if row["stratum"]=="published_positive" else a
    qc_z = z+pd.Timedelta(minutes=5) if row["stratum"]=="published_positive" else z
    total = int((qc_z-qc_a).total_seconds())+1
    valid = int(angles.loc[qc_a:qc_z].notna().sum())
    fraction = valid/total if total else 0.
    det = []
    center, stop = time(row["seed_start"]), time(row["seed_end"])
    for i, j, _, peak, cl, cr in v2.find_events(angles):
        x, y = b.index[i], b.index[j]
        if y < center-pd.Timedelta(minutes=15) or x > stop+pd.Timedelta(minutes=15):
            continue
        if not intersect(x, y+pd.Timedelta(seconds=1), a, z+pd.Timedelta(seconds=1)):
            continue
        det.append(dict(id=row["id"], method=method, start=x.isoformat(),
                        end=y.isoformat(), peak_angle=round(float(peak), 2),
                        left_censored=bool(cl), right_censored=bool(cr)))
    audit = dict(id=row["id"], method=method, status="ok" if fraction>=.8 else "low_coverage",
                 valid_fraction=round(fraction, 4), detections=len(det), error="")
    return audit, det, b


def raw_figure(row, b):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    mid = time(row["start"])+(time(row["end"])-time(row["start"]))/2
    a, z = mid-pd.Timedelta(minutes=15), mid+pd.Timedelta(minutes=15)
    d = b.loc[a-pd.Timedelta(minutes=10):z+pd.Timedelta(minutes=10)]
    fig, axs = plt.subplots(2,1,sharex=True,figsize=(12,6),constrained_layout=True)
    for comp in ("BR","BT","BN"):
        axs[0].plot(d.index, d[comp], label=comp, lw=.7)
    axs[1].plot(d.index, d.Bmag, color="black", lw=.8)
    for ax in axs:
        ax.axvspan(a, z, color="gray", alpha=.10)
        ax.grid(alpha=.2)
    axs[0].legend(ncol=3);axs[0].set_ylabel("RTN B (nT)")
    axs[1].set_ylabel("|B| (nT)");axs[1].set_xlabel("UTC")
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d %H:%M"))
    fig.suptitle(f"Independent raw FIELDS review {row['id']} (detector identities hidden)")
    destination=Path("figures/holdout_review")/f"{row['id']}.png"
    destination.parent.mkdir(exist_ok=True, parents=True)
    fig.savefig(destination,dpi=125);plt.close(fig)
    return str(destination)


def run(args):
    import reference_detector_v2 as v2
    manifest=OUT/"holdout_manifest.csv"
    settings=json.loads((OUT/"holdout_settings.json").read_text())
    if hashlib.sha256(manifest.read_bytes()).hexdigest()!=settings["manifest_sha256"]:
        raise SystemExit("Manifest modified after selection; abort!")
    selected=pd.read_csv(manifest,keep_default_na=False)
    audits=[];intervals=[];review=[]
    previous_path=OUT/"holdout_blind_review.csv"
    prior={}
    if previous_path.exists():
        old=pd.read_csv(previous_path, keep_default_na=False)
        prior={r["id"]:r for r in old.to_dict("records")}
    for i,row in enumerate(selected.to_dict("records"),1):
        b=None
        for method in METHODS:
            try:
                audit,events,magnetic=evaluate(row,method,v2)
                audits.append(audit)
                if audit["status"]=="ok": intervals.extend(events)
                if b is None: b=magnetic
                print(f"[{i}/{len(selected)}] {row['id']}/{method}: {audit['status']} ({len(events)} intervals)",flush=True)
            except Exception as exc:
                audits.append(dict(id=row["id"],method=method,status="error",
                                   valid_fraction=np.nan,detections=0,error=str(exc)[:200]))
                print(f"[{i}/{len(selected)}] {row['id']}/{method}: {exc}",flush=True)
        image=raw_figure(row,b) if b is not None and args.plots else ""
        mid=time(row["start"])+(time(row["end"])-time(row["start"]))/2
        previous=prior.get(row["id"],{})
        review.append(dict(id=row["id"],figure=image,
                           review_start=(mid-pd.Timedelta(minutes=15)).isoformat(),
                           review_end=(mid+pd.Timedelta(minutes=15)).isoformat(),
                           reviewer=previous.get("reviewer",""),
                           manual_classification=previous.get("manual_classification",""),
                           notes=previous.get("notes","")))
    a=pd.DataFrame(audits)
    det=pd.DataFrame(intervals,columns=["id","method","start","end","peak_angle","left_censored","right_censored"])
    a.to_csv(OUT/"holdout_audit.csv",index=False)
    det.to_csv(OUT/"holdout_detections.csv",index=False)
    pd.DataFrame(review).to_csv(previous_path,index=False)
    scores=[]
    for row in selected.to_dict("records"):
        statuses=a[a.id==row["id"]]
        if len(statuses)!=2 or not statuses.status.eq("ok").all():
            continue
        for method in METHODS:
            d=det[(det.id==row["id"]) & (det.method==method)]
            hit=bool(len(d))
            best=max((iou(time(row["start"]),time(row["end"]),time(x.start),time(x.end))
                      for x in d.itertuples(index=False)),default=0.)
            scores.append(dict(id=row["id"],stratum=row["stratum"],method=method,
                               candidate_count=len(d),any_candidate=hit,
                               best_iou=round(best,4) if row["stratum"]=="published_positive" else np.nan,
                               matched_iou025=best>=.25 if row["stratum"]=="published_positive" else None,
                               matched_iou050=best>=.5 if row["stratum"]=="published_positive" else None))
    score=pd.DataFrame(scores,columns=["id","stratum","method","candidate_count","any_candidate",
                                       "best_iou","matched_iou025","matched_iou050"])
    score.to_csv(OUT/"holdout_summary.csv",index=False)
    print("Paired successful windows:",score.id.nunique())
    for method in METHODS:
        for kind in ["published_positive","catalog_clear_unlabeled"]:
            g=score[(score.method==method)&(score.stratum==kind)]
            if g.empty:continue
            if kind=="published_positive":
                print(method,kind, "n=",len(g),"IoU>=0.25",int(g.matched_iou025.sum()),
                      "IoU>=0.50",int(g.matched_iou050.sum()))
            else:
                print(method,kind, "n=",len(g),"flagged",int(g.any_candidate.sum()),
                      "(NOT a measured false positive rate)")


def self_test(args=None):
    s=pd.Timestamp("2018-11-01")
    assert intersect(s,s+pd.Timedelta(minutes=5),s+pd.Timedelta(minutes=2),s+pd.Timedelta(minutes=8))
    assert iou(s,s+pd.Timedelta(minutes=10),s+pd.Timedelta(minutes=5),s+pd.Timedelta(minutes=15))==1/3
    refs=pd.DataFrame([dict(reference_id=f"E1-{i}", start=s+pd.Timedelta(hours=i*3+2),
                            end=s+pd.Timedelta(hours=i*3+2,minutes=2)) for i in range(20)])
    original=pd.DataFrame([dict(start="2018-11-01T03:00:00",end="2018-11-01T03:01:00")])
    coverage=pd.DataFrame([dict(day=x,valid_minutes=1440) for x in ["2018-11-01","2018-11-02","2018-11-03"]])
    one,_=sample(refs,original,coverage,4,2,11)
    two,_=sample(refs,original,coverage,4,2,11)
    assert len(one)==6 and one.equals(two)
    occupied=[span_from_seed(x) for x in original.to_dict("records")]
    for row in one.to_dict("records"):
        a,b=span_from_seed(dict(start=row["seed_start"],end=row["seed_end"]))
        assert not collides(a,b,occupied)
        occupied.append((a,b))
    print("PASS: deterministic holdout, no exposure reuse, interval IoU and control masking")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest="command",required=True)
    setup=sub.add_parser("prepare")
    setup.add_argument("--positives",type=int,default=10)
    setup.add_argument("--unlabeled",type=int,default=10)
    setup.add_argument("--seed",type=int,default=41801)
    execute=sub.add_parser("run")
    execute.add_argument("--plots",action="store_true")
    sub.add_parser("self-test")
    args=parser.parse_args()
    {"prepare":prepare,"run":run,"self-test":self_test}[args.command](args)
