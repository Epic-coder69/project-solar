"""Inspect NASA Parker Solar Probe SWEAP CDF variables and native sampling.

Run from repository root (see --help):
  python code/SWEAP/sweap_data_reader.py --list-products
  python code/SWEAP/sweap_data_reader.py --download --date 2018-11-06
  python code/SWEAP/sweap_data_reader.py --file data/sweap/spc-l3i/FILE.cdf --start 2018-11-06T12:00:00 --minutes 8
  python code/SWEAP/sweap_data_reader.py --self-test

SPC Level 3 contains derived moments and fits, NOT direct detector currents.
The Level 2 product contains collector-current and charge-flux spectra.
SWEAP time tags may be irregular and do not establish measurement exposure.
SPAN-I proton moments can be partial because the core falls outside its FOV.
NASA guide: https://spdf.gsfc.nasa.gov/pub/data/psp/sweap/sweap_data_user_guide.pdf
"""
import argparse
from datetime import date, timedelta
import json
from pathlib import Path
import re
import sys
from urllib.parse import urlparse

import numpy as np
import pandas as pd
import requests

CDAS = "https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/sp_phys/datasets"
PRODUCTS = {
    "spc-l3i": ("PSP_SWP_SPC_L3I", "psp_swp_spc_l3i", 0.2185,
                ["np_moment", "vp_moment_RTN", "wp_moment", "np1_fit", "vp1_fit_RTN", "wp1_fit"]),
    "spc-l2i": ("PSP_SWP_SPC_L2I", "psp_swp_spc_l2i", 0.2185,
                ["A_CURRENT", "B_CURRENT", "C_CURRENT", "D_CURRENT", "MV_LO", "MV_HI", "DIFF_CHARGE_FLUX_DENSITY"]),
    "spi-proton-l3": ("PSP_SWP_SPI_SF00_L3_MOM", "psp_swp_spi_sf00_l3_mom", None, []),
    "spe-pad-l3": ("PSP_SWP_SPE_SF0_L3_PAD", "psp_swp_spe_sf0_l3_pad", 13.981, []),
}


def utc(value):
    t = pd.Timestamp(value)
    return t.tz_localize("UTC") if t.tzinfo is None else t.tz_convert("UTC")


def cadence(times):
    """Actual adjacent EPOCH differences; NOT physical integration time."""
    a = pd.DatetimeIndex(times).dropna().sort_values()
    # Pandas 3 can preserve datetime64[s]/[ms]/[us] resolution. .asi8
    # reports values in that unit, not necessarily nanoseconds. Normalize
    # explicitly before converting adjacent EPOCH differences to seconds.
    diffs = np.diff(a.as_unit("ns").asi8) / 1e9
    positive = diffs[diffs > 0]
    return {"records": len(a), "duplicate_epochs": int(np.sum(diffs == 0)),
            "median_s": float(np.median(positive)) if len(positive) else None,
            "p10_s": float(np.percentile(positive, 10)) if len(positive) else None,
            "p90_s": float(np.percentile(positive, 90)) if len(positive) else None,
            "largest_gap_s": float(np.max(positive)) if len(positive) else None}


def valid(values, attributes):
    """Use variable-specific CDF fill and valid bounds; never replace bad data by 0."""
    a = np.asarray(values)
    if not np.issubdtype(a.dtype, np.number):
        return np.zeros(a.shape, dtype=bool)
    good = np.isfinite(a) & (np.abs(a) < 1e29)
    for key, op in (("FILLVAL", "fill"), ("VALIDMIN", "min"), ("VALIDMAX", "max")):
        if key in attributes:
            try:
                b = np.asarray(attributes[key], dtype=float)
                good &= (a != b if op == "fill" else a >= b if op == "min" else a <= b)
            except (TypeError, ValueError):
                pass
    return good


def discover(product, day, folder):
    """Get the most recent daily original NASA CDF, caching under git-ignored data/."""
    dataset, stem, _, _ = PRODUCTS[product]
    folder = folder / product
    folder.mkdir(parents=True, exist_ok=True)
    stamp = day.strftime("%Y%m%d")
    matcher = re.compile(rf"^{stem}_{stamp}_v(\d+)\.cdf$", re.I)
    candidates = [(int(m.group(1)), path) for path in folder.glob("*.cdf")
                  if (m := matcher.fullmatch(path.name))]
    if candidates:
        return max(candidates)[1]
    next_day = day + timedelta(days=1)
    endpoint = (f"{CDAS}/{dataset}/orig_data/"
                f"{day:%Y%m%d}T000000Z,{next_day:%Y%m%d}T000000Z/")
    response = requests.get(endpoint, headers={"Accept": "application/json"}, timeout=(20, 90))
    response.raise_for_status()
    descriptions = response.json().get("FileDescription", [])
    if isinstance(descriptions, dict):
        descriptions = [descriptions]
    entries = []
    for entry in descriptions:
        remote = entry.get("Name", "")
        name = Path(urlparse(remote).path).name
        if (match := matcher.fullmatch(name)) and remote.startswith("https://"):
            entries.append((int(match.group(1)), name, remote))
    if not entries:
        raise FileNotFoundError(f"No NASA {dataset} file for {day}; check mission coverage: {endpoint}")
    _, name, remote = max(entries)
    output = folder / name
    partial = folder / (name + ".part")
    print("NASA CDF:", remote, flush=True)
    try:
        with requests.get(remote, stream=True, timeout=(30, 180)) as stream:
            stream.raise_for_status()
            expected = stream.headers.get("Content-Length")
            count = 0
            with partial.open("wb") as target:
                for chunk in stream.iter_content(1024 * 1024):
                    if chunk:
                        target.write(chunk)
                        count += len(chunk)
        if count < 1024 or (expected and expected.isdecimal() and count != int(expected)):
            raise IOError(f"Incomplete or unexpectedly small CDF: {name}")
        partial.replace(output)
    except BaseException:
        partial.unlink(missing_ok=True)
        raise
    return output


def analyze(path, product, begin, minutes):
    try:
        import cdflib
    except ImportError as exc:
        raise RuntimeError("Install: python -m pip install -r code/FIELDS/external_libraries_requirements.txt") from exc
    cdf = cdflib.CDF(str(path))
    obj = cdf.cdf_info()
    names = list(obj.zVariables) + list(obj.rVariables)
    lookup = {name.lower(): name for name in names}
    if "epoch" not in lookup:
        raise ValueError(f"Cannot locate EPOCH in {names}")
    epoch = lookup["epoch"]
    times = pd.DatetimeIndex(pd.to_datetime(cdflib.cdfepoch.to_datetime(cdf.varget(epoch)), utc=True))
    good_time = times[~times.isna()]
    if good_time.empty:
        raise ValueError("CDF has no valid EPOCH records")
    start = utc(begin) if begin else good_time.min()
    stop = start + pd.Timedelta(minutes=minutes)
    indices = np.flatnonzero((times >= start) & (times < stop))
    frame = pd.DataFrame(index=times[indices]); frame.index.name = "UTC"
    report = {"source_file": str(path), "dataset": PRODUCTS[product][0],
              "product": product, "window_start_utc": start.isoformat(),
              "window_end_exclusive_utc": stop.isoformat(),
              "window_seconds": minutes * 60, "all_file_records": len(times),
              "all_file_cadence": cadence(times),
              "window_cadence": cadence(times[indices]),
              "nominal_catalog_cadence_seconds": PRODUCTS[product][2],
              "quality": {}, "variables": [],
              "notes": ["Epoch spacing is NOT an integration/exposure time or a proof of continuous coverage.",
                        "These are native samples, without time resampling or gap filling.",
                        "For science analysis, also inspect DQF flags, uncertainties, and FOV caveats."]}
    if product == "spi-proton-l3":
        report["notes"].append("SPAN-I proton moments are partial and can miss the core solar wind.")
    if product == "spc-l3i":
        report["notes"].append("SPC thermal speed principally describes the radial/reduced distribution; not full 3D temperature.")
    if product == "spc-l2i":
        report["notes"].append("SPC EPOCH is the start of a spectrum; individual steps have MEASUREMENT_TIME offsets.")
    q_name = lookup.get("general_flag")
    if q_name:
        flags = np.asarray(cdf.varget(q_name))
        flag_attrs = cdf.varattsget(q_name)
        good = valid(flags, flag_attrs) & (flags == 0)
        flag = flags[indices]
        frame["general_flag"] = flag
        report["quality"] = {"flag_variable": q_name,
            "good_zero": int(np.sum(good[indices])), "not_good_or_unknown": int(np.sum(~good[indices]))}
    else:
        good = None
        report["quality"] = {"warning": "No general_flag available. Cannot confirm good-quality samples."}
    wanted = PRODUCTS[product][3]
    if not wanted:
        wanted = [x for x in names if any(s in x.lower() for s in ("dens", "temp", "vel", "mom"))][:15]
    for requested in wanted:
        name = lookup.get(requested.lower())
        if name is None:
            continue
        attrs = cdf.varattsget(name)
        # Skip huge 3D distribution arrays in a lightweight discovery run.
        shape = list(getattr(cdf.varinq(name), "Dim_Sizes", []))
        detail = {"name": name, "units": str(attrs.get("UNITS", attrs.get("UNIT_PTR", ""))),
                  "description": str(attrs.get("CATDESC", ""))[:200], "dimensions": shape}
        if len(shape) > 1 or (shape and int(np.prod(shape)) > 128):
            detail["skipped"] = "Large multidimensional variable"
        else:
            values = np.asarray(cdf.varget(name))
            if values.ndim not in (1, 2) or len(values) != len(times):
                detail["skipped"] = "Not aligned to EPOCH"
            elif not np.issubdtype(values.dtype, np.number):
                detail["skipped"] = "Not numeric"
            else:
                selection = values[indices]
                mask = valid(selection, attrs)
                row_good = np.all(mask, axis=1) if selection.ndim == 2 else mask
                detail["valid_rows"] = int(np.sum(row_good))
                detail["valid_and_global_flag_zero"] = int(np.sum(row_good & good[indices])) if good is not None else None
                if selection.ndim == 1:
                    frame[name] = np.where(mask, selection, np.nan)
                else:
                    labels = ("R", "T", "N") if selection.shape[1] == 3 and "rtn" in name.lower() else tuple(map(str, range(selection.shape[1])))
                    for i, label in enumerate(labels):
                        frame[f"{name}_{label}"] = np.where(mask[:, i], selection[:, i], np.nan)
        report["variables"].append(detail)
    report["available_cdf_variables"] = names
    return report, frame


def test():
    data = np.array([5.0, -1e31, 6.0])
    assert valid(data, {"FILLVAL": -1e31}).tolist() == [True, False, True]
    times = pd.date_range("2018-11-06T00:00:00Z", periods=9, freq="1s").delete([3, 4])
    # Test mixed timestamp storage units: pandas 3 may use datetime64[s]
    # while CDF epochs generally arrive with micro/nanosecond resolution.
    for unit in ("s", "ms", "us", "ns"):
        got = cadence(times.as_unit(unit))
        assert got["largest_gap_s"] == 3.0, (unit, got)
        assert got["median_s"] == 1.0, (unit, got)
        assert got["records"] == 7, (unit, got)
    nanos = pd.DatetimeIndex(["2018-11-06T00:00:00.000000001Z",
                              "2018-11-06T00:00:00.873820033Z"])
    assert abs(cadence(nanos)["median_s"] - 0.873820032) < 1e-8
    assert utc("2018-11-06T00:00:00").tzinfo is not None
    print("PASS: offline fill-value, cadence, gap, and UTC checks")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--list-products", action="store_true")
    p.add_argument("--product", default="spc-l3i", choices=PRODUCTS)
    p.add_argument("--date", help="UTC YYYY-MM-DD (required with --download)")
    p.add_argument("--download", action="store_true")
    p.add_argument("--file", type=Path, help="Existing local CDF instead of download")
    p.add_argument("--start", help="Window UTC start, default earliest file EPOCH")
    p.add_argument("--minutes", type=float, default=8)
    p.add_argument("--cache-dir", type=Path, default=Path("data/sweap"))
    p.add_argument("--csv", type=Path, help="Optional CSV with native window measurements")
    p.add_argument("--json", type=Path, help="Optional JSON diagnostic report")
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    if args.self_test:
        test(); return
    if args.list_products:
        for key, (dataset, _, cadence_s, _) in PRODUCTS.items():
            print(f"{key:16} {dataset:32} advertised cadence {cadence_s} s")
        return
    if not np.isfinite(args.minutes) or args.minutes <= 0:
        p.error("--minutes must be positive")
    if bool(args.file) == bool(args.download):
        p.error("Choose exactly one: --file CDF or --download --date YYYY-MM-DD")
    if args.download and not args.date:
        p.error("--download requires --date")
    path = args.file or discover(args.product, date.fromisoformat(args.date), args.cache_dir)
    if not path.is_file():
        p.error(f"File does not exist: {path}")
    report, data = analyze(path, args.product, args.start, args.minutes)
    print(f"\nDataset {report['dataset']}\nFile {path}")
    print(f"Window [{report['window_start_utc']}, {report['window_end_exclusive_utc']})")
    print("Full file EPOCH cadence:", report["all_file_cadence"])
    print("8-minute window EPOCH cadence:", report["window_cadence"])
    print("Quality:", report["quality"])
    for item in report["variables"]:
        print(" -", item)
    print("All NASA CDF variable names:", ", ".join(report["available_cdf_variables"]))
    print("\nSample native records (not resampled):")
    print(data.head(8).to_string())
    print("\nInterpretation cautions:")
    for note in report["notes"]:
        print(" -", note)
    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        data.to_csv(args.csv)
        print("Saved:", args.csv)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print("Saved:", args.json)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, requests.RequestException, RuntimeError) as e:
        sys.exit(f"SWEAP inspection failed: {e}")
