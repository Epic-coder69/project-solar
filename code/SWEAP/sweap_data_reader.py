"""Unified NASA Parker Solar Probe SWEAP CDF reader and family inspector.

Run from repository root (see --help):
  python code/SWEAP/sweap_data_reader.py --list-products
  python code/SWEAP/sweap_data_reader.py --survey --date 2018-11-06 --catalog-only
  python code/SWEAP/sweap_data_reader.py --survey --date 2018-11-06 --subset-on-limit
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



# ---- NASA SWEAP eight-family CDF catalog and bounded preview utilities ----

# Family, real CDAWeb dataset ID, exact NASA filename stem.
CATALOG = {
    "spc-l2i": ("cup-charge-flux", "PSP_SWP_SPC_L2I", "psp_swp_spc_l2i"),
    "spc-l3i": ("proton-moments", "PSP_SWP_SPC_L3I", "psp_swp_spc_l3i"),
    "spi-proton-l2": ("proton-3d", "PSP_SWP_SPI_SF00_L2_8DX32EX8A", "psp_swp_spi_sf00_l2_8dx32ex8a"),
    "spi-alpha-l2": ("alpha-3d", "PSP_SWP_SPI_SF01_L2_8DX32EX8A", "psp_swp_spi_sf01_l2_8dx32ex8a"),
    "spi-proton-l3": ("ion-moments", "PSP_SWP_SPI_SF00_L3_MOM", "psp_swp_spi_sf00_l3_mom"),
    "spi-alpha-l3": ("ion-moments", "PSP_SWP_SPI_SF0A_L3_MOM", "psp_swp_spi_sf0a_l3_mom"),
    "spa-electron-3d": ("electron-3d", "PSP_SWP_SPA_SF0_L2_16AX8DX32E", "psp_swp_spa_sf0_l2_16ax8dx32e"),
    "spb-electron-3d": ("electron-3d", "PSP_SWP_SPB_SF0_L2_16AX8DX32E", "psp_swp_spb_sf0_l2_16ax8dx32e"),
    "spa-electron-spectrum": ("electron-energy", "PSP_SWP_SPA_SF1_L2_32E", "psp_swp_spa_sf1_l2_32e"),
    "spb-electron-spectrum": ("electron-energy", "PSP_SWP_SPB_SF1_L2_32E", "psp_swp_spb_sf1_l2_32e"),
    "spe-pad-l3": ("electron-pitch-angle", "PSP_SWP_SPE_SF0_L3_PAD", "psp_swp_spe_sf0_l3_pad"),
    "spa-pad-l3": ("electron-pitch-angle", "PSP_SWP_SPA_SF0_L3_PAD", "psp_swp_spa_sf0_l3_pad"),
    "spb-pad-l3": ("electron-pitch-angle", "PSP_SWP_SPB_SF0_L3_PAD", "psp_swp_spb_sf0_l3_pad"),
}

CAVEATS = {
    "cup-charge-flux": "SPC L2 is calibrated charge-flux/current spectra, not raw telemetry. Inspect per-step timing.",
    "proton-moments": "SPC L3 contains derived moments/fits; adjacent timestamps do not imply full exposure.",
    "proton-3d": "SPAN-I field of view can exclude the solar-wind proton core.",
    "alpha-3d": "SF01 alpha flux is PROTON-CONTAMINATED; do not treat it as corrected alpha flux.",
    "ion-moments": "SPAN-I proton/alpha moments are partial due to field-of-view limitations.",
    "electron-3d": "SPAN-E 2018 FOV obstructions, low-energy secondary electrons and invalid pixels matter.",
    "electron-energy": "SF1 spectra sum over angles; high-deflection bias makes them unsuitable for rigid statistics.",
    "electron-pitch-angle": "PAD depends on FIELDS magnetic direction; check electron quality/coverage.",
}


def file_for_day(entries, stem, day):
    """Choose the newest exact-day version; never silently use a nearby day."""
    pattern = re.compile(rf"^{re.escape(stem)}_{day:%Y%m%d}_v(\d+)\.cdf$", re.I)
    found = []
    for entry in entries:
        remote = entry.get("Name", "")
        name = Path(urlparse(remote).path).name
        match = pattern.fullmatch(name)
        if match and remote.startswith("https://"):
            found.append((int(match.group(1)), name, remote))
    return max(found, default=None)


def available_file(product, day, cache):
    _, dataset, stem = CATALOG[product]
    folder = cache / product
    folder.mkdir(parents=True, exist_ok=True)
    cached = []
    for path in folder.glob("*.cdf"):
        match = file_for_day([{"Name": "https://cached/" + path.name}], stem, day)
        if match:
            cached.append((match[0], str(path)))
    if cached:
        return {"status": "cached", "path": max(cached)[1]}
    tomorrow = day + timedelta(days=1)
    endpoint = (f"{CDAS}/{dataset}/orig_data/"
                f"{day:%Y%m%d}T000000Z,{tomorrow:%Y%m%d}T000000Z/")
    response = requests.get(endpoint, timeout=45, headers={"Accept": "application/json"})
    response.raise_for_status()
    entries = response.json().get("FileDescription", [])
    if isinstance(entries, dict):
        entries = [entries]
    chosen = file_for_day(entries, stem, day)
    if not chosen:
        return {"status": "not_available_on_date", "catalog_url": endpoint}
    version, filename, remote = chosen
    return {"status": "found", "path": str(folder / filename),
            "url": remote, "version": version}


def bounded_download(url, path, max_bytes, remaining_bytes):
    """Never exceed file/total size caps, including unknown-length responses."""
    if min(max_bytes, remaining_bytes) <= 0:
        return "skipped_size_budget", 0
    path = Path(path)
    partial = path.with_suffix(path.suffix + ".part")
    try:
        with requests.get(url, stream=True, timeout=(30, 120)) as response:
            response.raise_for_status()
            expected = response.headers.get("Content-Length")
            limit = min(max_bytes, remaining_bytes)
            if expected and expected.isdigit() and int(expected) > limit:
                return "skipped_size_budget", 0
            count = 0
            with partial.open("wb") as output:
                for chunk in response.iter_content(1024 * 1024):
                    if chunk:
                        count += len(chunk)
                        if count > limit:
                            return "skipped_size_budget", 0
                        output.write(chunk)
        if count < 1024 or (expected and expected.isdigit() and count != int(expected)):
            return "incomplete_download", 0
        partial.replace(path)
        return "downloaded", count
    finally:
        partial.unlink(missing_ok=True)


def subset_for_window(product, day, cache, start, minutes):
    """Request a NASA-generated, native-resolution CDF for just the UTC window.

    This is a server-side time subset (NOT averaging/resampling). Never mix
    these files with the official original daily CDFs. The CDAWeb REST
    service can return errors or decline to subset particular large products.
    """
    dataset = CATALOG[product][1]
    begin = pd.Timestamp(start or day.isoformat())
    begin = begin.tz_localize("UTC") if begin.tzinfo is None else begin.tz_convert("UTC")
    stop = begin + pd.Timedelta(minutes=minutes)
    # REST endpoint has whole-second ISO 8601 timestamps; round outwards.
    request_begin = begin.floor("s")
    request_stop = stop.ceil("s")
    time_label = lambda t: t.strftime("%Y%m%dT%H%M%SZ")
    lo, hi = time_label(request_begin), time_label(request_stop)
    path = cache / "_subsets" / product / f"{product}_{lo}_{hi}.cdf"
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file():
        return {"status": "cached", "path": str(path), "mode": "time_subset"}
    endpoint = f"{CDAS}/{dataset}/data/{lo},{hi}/ALL-VARIABLES?format=cdf"
    response = requests.get(endpoint, timeout=(20, 120),
                            headers={"Accept": "application/json"})
    response.raise_for_status()
    payload = response.json()
    descriptions = payload.get("FileDescription", [])
    if isinstance(descriptions, dict):
        descriptions = [descriptions]
    # Do not download an arbitrary webpage or an ambiguous multi-file result.
    urls = [entry.get("Name", "") for entry in descriptions
            if entry.get("Name", "").startswith("https://")
            and urlparse(entry.get("Name", "")).path.lower().endswith(".cdf")]
    if len(urls) != 1:
        message = payload.get("Status", payload.get("CdasStatus", {}))
        return {"status": "subset_unavailable", "catalog_url": endpoint,
                "reason": f"Expected one CDF link, got {len(urls)}; status={str(message)[:200]}"}
    return {"status": "found", "path": str(path), "url": urls[0],
            "mode": "time_subset", "request_url": endpoint}


def _epoch(cdf):
    import cdflib
    info = cdf.cdf_info()
    names = list(info.zVariables) + list(info.rVariables)
    key = next((v for v in names if v.lower() == "epoch"), None)
    if key is None:
        key = next((v for v in names if v.lower().startswith("epoch")), None)
    if key is None:
        return None, None
    values = cdf.varget(key)
    times = pd.DatetimeIndex(pd.to_datetime(cdflib.cdfepoch.to_datetime(values), utc=True))
    return key, times


def _peek(cdf, name, row, max_examples=6):
    """Read only one time record of potentially enormous 4D CDF arrays."""
    metadata = cdf.varinq(name)
    if int(metadata.Last_Rec) < row or int(metadata.Last_Rec) < 0:
        return None
    sample = np.asarray(cdf.varget(name, startrec=row, endrec=row))
    first = sample.ravel()[:max_examples]
    if np.issubdtype(sample.dtype, np.number):
        values = [float(x) if np.isfinite(x) else None for x in first]
    else:
        values = [str(x)[:80] for x in first]
    return {"record": row, "returned_shape": list(sample.shape), "values": values}


def inspect(path, window_start, minutes, variable_limit=20):
    """Report axes, units, variable names, representative records and spacing."""
    import cdflib
    cdf = cdflib.CDF(str(path))
    info = cdf.cdf_info()
    names = list(info.zVariables) + list(info.rVariables)
    epoch, times = _epoch(cdf)
    result = {"all_variable_names": names, "epoch_variable": epoch,
              "variable_count": len(names), "variables": [],
              "quality_flag_names": [x for x in names if "flag" in x.lower() or "dqf" in x.lower()]}
    row = 0
    window_has_records = False
    if times is not None and len(times):
        start = pd.Timestamp(window_start)
        start = start.tz_localize("UTC") if start.tzinfo is None else start.tz_convert("UTC")
        end = start + pd.Timedelta(minutes=minutes)
        ids = np.flatnonzero((times >= start) & (times < end))
        if len(ids):
            row = int(ids[0])
            window_has_records = True
        delta = np.diff(times.as_unit("ns").asi8) / 1e9
        delta = delta[delta > 0]
        within = np.diff(times[ids].as_unit("ns").asi8) / 1e9 if len(ids) > 1 else np.array([])
        result["timing"] = {"file_records": len(times), "window_records": len(ids),
            "start_utc": start.isoformat(), "window_minutes": minutes,
            "file_median_epoch_spacing_s": float(np.median(delta)) if len(delta) else None,
            "window_median_epoch_spacing_s": float(np.median(within)) if len(within) else None,
            "window_max_epoch_gap_s": float(np.max(within)) if len(within) else None}
    highlights = ("flux", "current", "dens", "moment", "temp", "vel", "fit",
                  "pad", "energy", "quality", "dqf", "theta", "phi", "epoch")
    ordered = sorted(names, key=lambda n: (not any(x in n.lower() for x in highlights), n.lower()))
    for name in ordered[:variable_limit]:
        try:
            v = cdf.varinq(name)
            attrs = cdf.varattsget(name)
            entry = {"name": name, "dimensions_per_record": list(v.Dim_Sizes),
                     "last_record": int(v.Last_Rec), "units": str(attrs.get("UNITS", "")),
                     "description": str(attrs.get("CATDESC", ""))[:160],
                     "depend_0": str(attrs.get("DEPEND_0", "")),
                     "var_type": str(attrs.get("VAR_TYPE", ""))}
            # Never show an outside-window record as an example of this window.
            if name != epoch and window_has_records and int(v.Last_Rec) >= row:
                entry["example"] = _peek(cdf, name, row)
            result["variables"].append(entry)
        except (ValueError, OSError, TypeError, AttributeError) as exc:
            result["variables"].append({"name": name, "read_error": str(exc)})
    return result


def survey(day, cache, output, *, minutes=8, start=None, catalog_only=False,
           max_file_mib=50, max_total_mib=250, families=None, subset_on_limit=False):
    """Survey every family; report missing/oversized/error products separately."""
    selected = [(name, cfg) for name, cfg in CATALOG.items()
                if not families or cfg[0] in families]
    results = {"date": day.isoformat(), "window_minutes": minutes,
               "mode": "catalog_only" if catalog_only else "bounded_cdf_inspection",
               "subset_on_limit": subset_on_limit,
               "families": sorted({cfg[0] for _, cfg in selected}), "products": [],
               "caveats": CAVEATS, "limitations": [
                   "Found in NASA catalog is not the same as observed valid samples.",
                   "CDF EPOCH spacing is not physical integration/exposure duration.",
                   "All examples are native record snippets, not averaged or resampled.",
                   "Quality bit meanings, fill values, FOV, and energy axes require further study."]}
    remaining = int(max_total_mib * 1048576)
    for name, (family, dataset, _) in selected:
        item = {"product": name, "family": family, "dataset": dataset,
                "interpretation_caveat": CAVEATS[family]}
        try:
            found = available_file(name, day, cache)
            item.update({"status": found["status"], "file": found.get("path"),
                         "source_url": found.get("url")})
            if item["status"] == "found" and not catalog_only:
                state, size = bounded_download(found["url"], found["path"],
                                               int(max_file_mib * 1048576), remaining)
                remaining -= size
                item["status"] = state
                item["bytes_downloaded"] = size
            if item["status"] == "skipped_size_budget" and subset_on_limit and not catalog_only:
                subset = subset_for_window(name, day, cache, start, minutes)
                item["original_daily_status"] = "skipped_size_budget"
                item["subset_status"] = subset["status"]
                item["subset_request_url"] = subset.get("request_url", subset.get("catalog_url"))
                if subset["status"] in ("found", "cached"):
                    found = subset
                    item["file"] = found["path"]
                    item["source_url"] = found.get("url", "local_cache")
                    item["mode"] = "time_subset"
                    if subset["status"] == "found":
                        state, size = bounded_download(found["url"], found["path"],
                                                       int(max_file_mib * 1048576), remaining)
                        remaining -= size
                        item["bytes_downloaded"] = size
                        item["status"] = state
                    else:
                        item["status"] = "cached"
                else:
                    item["status"] = subset["status"]
                    item["reason"] = subset.get("reason")
            if item["status"] in ("cached", "downloaded") and not catalog_only:
                item["cdf"] = inspect(found["path"], start or day.isoformat(), minutes)
                item["status"] = "inspected"
        except (requests.RequestException, ValueError, OSError, ImportError, KeyError) as exc:
            item["status"] = "error"
            item["error"] = f"{type(exc).__name__}: {exc}"
        results["products"].append(item)
        rows = item.get("cdf", {}).get("timing", {}).get("window_records", "-")
        print(f"[{family:21}] {name:25} {item['status']:24} window_records={rows}", flush=True)
    results["family_status"] = {
        family: {"inspected_products": [p["product"] for p in results["products"]
                                        if p["family"] == family and p["status"] == "inspected"],
                 "product_statuses": {p["product"]: p["status"] for p in results["products"]
                                      if p["family"] == family}}
        for family in results["families"]
    }
    results["summary"] = {
        "inspected": sum(p["status"] == "inspected" for p in results["products"]),
        "not_available": sum(p["status"] == "not_available_on_date" for p in results["products"]),
        "size_limited": sum(p["status"] == "skipped_size_budget" for p in results["products"]),
        "subset_unavailable": sum(p["status"] == "subset_unavailable" for p in results["products"]),
        "errors": sum(p["status"] == "error" for p in results["products"])}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(results, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"[SWEAP SURVEY] {results['summary']}\n[SAVED] {output}")
    return results


def self_test():
    day = date(2018, 11, 6)
    entries = [{"Name": "https://nasa/files/psp_swp_spc_l3i_20181106_v01.cdf"},
               {"Name": "https://nasa/files/psp_swp_spc_l3i_20181106_v22.cdf"},
               {"Name": "https://nasa/files/psp_swp_spc_l3i_20181107_v99.cdf"}]
    chosen = file_for_day(entries, "psp_swp_spc_l3i", day)
    assert chosen[0] == 22 and chosen[1].endswith("v22.cdf")
    assert len(CATALOG) == 13 and len({v[0] for v in CATALOG.values()}) == 8
    assert set(CAVEATS) == {v[0] for v in CATALOG.values()}
    # A complete 3D flux day could be enormous. Confirm that previewing
    # it only asks cdflib for ONE native time record and preserves axes.
    from types import SimpleNamespace
    class FakeCDF:
        def varinq(self, name):
            return SimpleNamespace(Last_Rec=100, Dim_Sizes=[8, 32, 8])
        def varget(self, name, *, startrec, endrec):
            assert startrec == endrec == 4
            return np.arange(2048).reshape((1, 8, 32, 8))
    sample = _peek(FakeCDF(), "synthetic_proton_3d_flux", 4)
    assert sample["returned_shape"] == [1, 8, 32, 8]
    assert sample["values"] == [0., 1., 2., 3., 4., 5.]
    # Test the size cap without any real network traffic.
    from unittest.mock import patch
    from tempfile import TemporaryDirectory
    class FakeResponse:
        headers = {"Content-Length": "99999999"}
        def __enter__(self):
            return self
        def __exit__(self, *args):
            return False
        def raise_for_status(self):
            pass
        def iter_content(self, _):
            raise AssertionError("oversized download must not begin")
    with TemporaryDirectory() as folder, patch.object(requests, "get", return_value=FakeResponse()):
        state, count = bounded_download("https://example.invalid/test.cdf",
                                        Path(folder) / "test.cdf", 4096, 4096)
        assert state == "skipped_size_budget" and count == 0
        assert not list(Path(folder).glob("*.part"))
    # NASA's CDAS REST can return a URL for a CDF containing only the
    # selected time interval. Verify its UTC bounds and separate cache path.
    class FakeSubsetResponse:
        def raise_for_status(self):
            pass
        def json(self):
            return {"FileDescription": [{
                "Name": "https://cdaweb.gsfc.nasa.gov/sp_phys/tmp/preview.cdf"}]}
    with TemporaryDirectory() as folder, patch.object(requests, "get",
                                                     return_value=FakeSubsetResponse()) as get:
        subset = subset_for_window("spi-proton-l2", day, Path(folder),
                                   "2018-11-06T00:00:00Z", 8)
        assert subset["status"] == "found"
        assert subset["mode"] == "time_subset"
        assert "_subsets" in subset["path"]
        assert "20181106T000000Z_20181106T000800Z" in subset["path"]
        assert "ALL-VARIABLES?format=cdf" in get.call_args.args[0]
    print("PASS: SWEAP eight-family catalog, 13 product IDs, exact-day selection, "
          "3D one-record sampling, strict download budgets, and NASA time-subset request")



# --- FIELDS-driven native SWEAP event extraction --------------------------

def load_event_intervals(path, max_events=0, before_seconds=10, after_seconds=10):
    """Read exact FIELDS half-open event intervals, with optional UTC flanks."""
    import hashlib
    frame = pd.read_csv(path, dtype={"event_id": str})
    required = {"event_start_utc", "event_end_utc"}
    if not required.issubset(frame):
        raise ValueError("FIELDS CSV missing " + ", ".join(sorted(required - set(frame.columns))))
    if frame.empty:
        raise ValueError("FIELDS event CSV contains no events")
    if max_events:
        frame = frame.iloc[:max_events]
    result = []
    for i, row in frame.iterrows():
        start = utc(row["event_start_utc"])
        if "event_end_exclusive_utc" in frame and pd.notna(row["event_end_exclusive_utc"]):
            end = utc(row["event_end_exclusive_utc"])
        else:
            # FIELDS event_end_utc marks the inclusive one-second magnetic bin.
            end = utc(row["event_end_utc"]) + pd.Timedelta(seconds=1)
        if end <= start:
            raise ValueError(f"FIELDS event row {i}: exclusive end must follow start")
        event_id = (str(row["event_id"]) if "event_id" in frame
                    and pd.notna(row["event_id"]) else f"row_{i+1}")
        begin = start - pd.Timedelta(seconds=before_seconds)
        stop = end + pd.Timedelta(seconds=after_seconds)
        # No user-provided event ID appears in output paths (path traversal safe).
        key = hashlib.sha256(f"{event_id}|{start.isoformat()}|{end.isoformat()}".encode()).hexdigest()[:16]
        result.append({"event_id": event_id, "event_key": key,
                       "event_start_utc": start.isoformat(),
                       "event_end_exclusive_utc": end.isoformat(),
                       "window_start_utc": begin.isoformat(),
                       "window_end_exclusive_utc": stop.isoformat()})
    return result


def split_utc_days(begin, stop):
    """Split [begin,stop) at UTC midnight with no duplicated boundaries."""
    a, b = utc(begin), utc(stop)
    while a < b:
        next_stop = min(b, a.floor("D") + pd.Timedelta(days=1))
        yield a.date(), a, next_stop
        a = next_stop


def native_event_slice(cdf, begin, stop, max_variable_mib=128):
    """Extract native CDF records from all variables with provable time axes.

    Preserves multidimensional arrays, detector flags and fill values.
    Never interpolates, rebins, averages or substitutes missing records.
    """
    import cdflib
    info = cdf.cdf_info()
    names = list(info.zVariables) + list(info.rVariables)
    lookup = {n.lower(): n for n in names}
    timestamps = {}
    for name in names:
        if not name.lower().startswith("epoch"):
            continue
        try:
            raw = cdf.varget(name)
            times = cdflib.cdfepoch.to_datetime(raw)
            timestamps[name] = pd.DatetimeIndex(pd.to_datetime(times, utc=True))
        except (ValueError, TypeError, OverflowError):
            continue
    if not timestamps:
        raise ValueError("CDF has no readable UTC Epoch/Epoch* axis")
    result = {"source_variable_count": len(names), "variables": {}, "errors": {},
              "time_axes": {}, "unassigned_record_varying": [],
              "full_native_arrays": True}
    selections = {}
    left, right = utc(begin), utc(stop)
    for name, times in timestamps.items():
        ids = np.flatnonzero((times >= left) & (times < right))
        selections[name] = ids
        spacing = (np.diff(times[ids].as_unit("ns").asi8) / 1e9
                   if len(ids) > 1 else np.array([]))
        result["time_axes"][name] = {
            "records_in_window": int(len(ids)), "file_records": int(len(times)),
            "median_epoch_step_s": float(np.median(spacing)) if len(spacing) else None,
            "max_epoch_gap_s": float(np.max(spacing)) if len(spacing) else None,
            "first_utc": times[ids[0]].isoformat() if len(ids) else None,
            "last_utc": times[ids[-1]].isoformat() if len(ids) else None}
    arrays = {}
    max_bytes = int(max_variable_mib * 1048576)
    for name in names:
        try:
            var = cdf.varinq(name)
            attrs = cdf.varattsget(name)
            varying = bool(var.Rec_Vary)
            nrecords = int(var.Last_Rec) + 1
            dep = str(attrs.get("DEPEND_0", "")).strip()
            if dep not in timestamps:
                dep = lookup.get(dep.lower(), dep)
            inferred = False
            if not varying:
                rowids = None
            else:
                if dep not in timestamps:
                    candidates = [ax for ax, times in timestamps.items()
                                  if len(times) == nrecords]
                    if len(candidates) == 1:
                        dep = candidates[0]
                        inferred = True
                    else:
                        result["unassigned_record_varying"].append(name)
                        result["variables"][name] = {
                            "status": "skipped_unresolved_time_axis",
                            "rec_vary": True, "depend_0": dep,
                            "records": nrecords}
                        continue
                rowids = selections[dep]
                if not len(rowids):
                    result["variables"][name] = {
                        "status": "zero_records_in_window",
                        "rec_vary": True, "depend_0": dep}
                    continue
                if int(rowids[-1]) >= nrecords:
                    result["variables"][name] = {
                        "status": "skipped_inconsistent_time_axis",
                        "rec_vary": True, "depend_0": dep,
                        "records": nrecords}
                    continue
            if varying:
                lo, hi = int(rowids[0]), int(rowids[-1])
                slab = np.asarray(cdf.varget(name, startrec=lo, endrec=hi))
                slab = slab[(rowids - lo).astype(int)]
            else:
                slab = np.asarray(cdf.varget(name))
            if slab.nbytes > max_bytes:
                result["variables"][name] = {
                    "status": "skipped_variable_size_limit",
                    "bytes": int(slab.nbytes), "max_bytes": max_bytes}
                continue
            if slab.dtype.kind == "O":
                # NPZ files must never require pickle. Plain string objects
                # are safe to convert to Unicode; others are explicitly omitted.
                if any(not isinstance(v, (str, bytes)) for v in slab.ravel()):
                    result["variables"][name] = {"status": "skipped_object_dtype"}
                    continue
                slab = slab.astype(str)
            arrays[name] = slab
            result["variables"][name] = {
                "status": "saved", "rec_vary": varying,
                "depend_0": dep if varying else None,
                "time_axis_inferred": inferred,
                "shape": list(slab.shape), "dtype": str(slab.dtype),
                "units": str(attrs.get("UNITS", "")),
                "catdesc": str(attrs.get("CATDESC", ""))[:200],
                "fillval": str(attrs.get("FILLVAL", "")),
                "bytes": int(slab.nbytes)}
        except (ValueError, OSError, TypeError, IndexError, AttributeError) as exc:
            result["variables"][name] = {
                "status": "read_error", "reason": f"{type(exc).__name__}: {exc}"}
            result["errors"][name] = f"{type(exc).__name__}: {exc}"
    result["saved_variables"] = sum(v["status"] == "saved"
                                    for v in result["variables"].values())
    result["not_saved_variables"] = len(names) - result["saved_variables"]
    result["any_records_in_window"] = any(len(v) for v in selections.values())
    return result, arrays


def event_source(product, day, start, stop, cache, state, budget, file_mib):
    """Reuse a daily CDF when small, else request a NASA time subset."""
    key = (product, str(day))
    if key not in state:
        try:
            found = available_file(product, day, cache)
            if found["status"] == "found":
                status, used = bounded_download(found["url"], found["path"],
                                                int(file_mib * 1048576), budget[0])
                budget[0] -= used
                found["status"] = status
                if status == "skipped_size_budget":
                    found["reason"] = "daily_file_exceeds_budget"
            state[key] = found
        except (requests.RequestException, OSError, ValueError) as exc:
            state[key] = {"status": "error", "reason": str(exc)}
    source = state[key]
    if source["status"] in ("cached", "downloaded"):
        return source
    if source["status"] != "skipped_size_budget":
        return source
    if budget[0] <= 0:
        return {"status": "skipped_total_size_budget"}
    try:
        minutes = float((stop - start) / pd.Timedelta(minutes=1))
        subset = subset_for_window(product, day, cache, start.isoformat(), minutes)
        if subset["status"] == "found":
            status, used = bounded_download(subset["url"], subset["path"],
                                            int(file_mib * 1048576), budget[0])
            budget[0] -= used
            subset["status"] = status
        return subset
    except (requests.RequestException, OSError, ValueError) as exc:
        return {"status": "error", "reason": str(exc)}


def extract_events(csv_path, output_dir, cache, products, *, max_events=0,
                   before_seconds=10, after_seconds=10,
                   max_file_mib=50, max_total_mib=2048,
                   max_variable_mib=128, resume=True):
    """Write all available native measurements plus an explicit coverage manifest."""
    import cdflib
    events = load_event_intervals(csv_path, max_events, before_seconds, after_seconds)
    output_dir.mkdir(parents=True, exist_ok=True)
    budget = [int(max_total_mib * 1048576)]
    daily_state = {}
    rows = []
    print(f"[SWEAP EVENTS] {len(events)} FIELDS candidates; {len(products)} NASA products; "
          f"flanks {before_seconds}/{after_seconds}s", flush=True)
    for product in products:
        for event in events:
            for day, begin, end in split_utc_days(event["window_start_utc"],
                                                  event["window_end_exclusive_utc"]):
                path = (output_dir / product /
                        f"{event['event_key']}_{day:%Y%m%d}_{begin:%H%M%S%f}.npz")
                meta_path = path.with_suffix(".json")
                row = {**event, "product": product, "dataset": CATALOG[product][1],
                       "day": day.isoformat(),
                       "segment_start_utc": begin.isoformat(),
                       "segment_end_exclusive_utc": end.isoformat(),
                       "archive": str(path), "metadata": str(meta_path)}
                if resume and path.is_file() and meta_path.is_file():
                    try:
                        prior = json.loads(meta_path.read_text(encoding="utf-8"))
                        if prior.get("status") == "saved":
                            row.update({"status": "resumed",
                                        "saved_variables": prior.get("saved_variables", 0),
                                        "records_primary_epoch": prior.get("records_primary_epoch", 0)})
                            rows.append(row)
                            continue
                    except (OSError, ValueError):
                        pass
                src = event_source(product, day, begin, end, cache, daily_state,
                                   budget, max_file_mib)
                row["source_status"] = src["status"]
                row["source_file"] = src.get("path")
                row["source_url"] = src.get("url")
                row["source_kind"] = src.get("mode", "original_daily_cdf")
                if src["status"] not in ("cached", "downloaded"):
                    row.update(status="no_source", reason=src.get("reason", src["status"]))
                else:
                    try:
                        report, arrays = native_event_slice(
                            cdflib.CDF(str(src["path"])), begin, end, max_variable_mib)
                        row["saved_variables"] = report["saved_variables"]
                        primary = next((n for n in report["time_axes"] if n.lower() == "epoch"),
                                       next(iter(report["time_axes"])))
                        row["records_primary_epoch"] = report["time_axes"][primary]["records_in_window"]
                        row["status"] = ("saved" if report["any_records_in_window"]
                                         else "no_epoch_records")
                        row["not_saved_variables"] = report["not_saved_variables"]
                        row["read_errors"] = len(report["errors"])
                        row["unassigned_time_variables"] = report["unassigned_record_varying"]
                        if row["status"] == "saved":
                            path.parent.mkdir(parents=True, exist_ok=True)
                            tmp = path.with_suffix(".part.npz")
                            try:
                                np.savez_compressed(tmp, **arrays)
                                tmp.replace(path)
                            finally:
                                tmp.unlink(missing_ok=True)
                        row["variable_details"] = report
                    except (ValueError, OSError, TypeError, IndexError, AttributeError) as exc:
                        row.update(status="read_error",
                                   reason=f"{type(exc).__name__}: {exc}")
                meta_path.parent.mkdir(parents=True, exist_ok=True)
                meta_path.write_text(json.dumps(row, indent=2, default=str) + "\n",
                                     encoding="utf-8")
                rows.append({k: v for k, v in row.items() if k != "variable_details"})
        totals = pd.Series([r["status"] for r in rows if r["product"] == product]).value_counts().to_dict()
        print(f"[SWEAP PRODUCT] {product}: {totals}", flush=True)
    manifest = output_dir / "event_manifest.csv"
    pd.DataFrame(rows).to_csv(manifest, index=False)
    summary = {
        "events_requested": len(events), "products_requested": len(products),
        "segments": len(rows),
        "status_counts": pd.Series([r["status"] for r in rows]).value_counts().to_dict(),
        "bytes_downloaded": int(max_total_mib * 1048576) - budget[0],
        "event_source_csv": str(csv_path), "manifest_csv": str(manifest),
        "notes": [
            "Native CDF records and dimensions preserved without resampling.",
            "Epoch time within an event does not establish measurement exposure overlap.",
            "Unaligned variables and read errors are listed explicitly.",
            "Saved records may include fill values and flags; none are assumed science-valid."]}
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n",
                                               encoding="utf-8")
    print(f"[SWEAP EVENTS SAVED] {summary['status_counts']}\n[MANIFEST] {manifest}",
          flush=True)
    return summary



def event_self_test():
    """Offline regression test of exact FIELDS intervals and full 3D records."""
    import tempfile
    from types import SimpleNamespace, ModuleType
    from unittest.mock import patch

    class FakeCDF:
        def __init__(self):
            self.Epoch = np.arange(20, dtype=np.int64)
            self.Epoch_ALT = np.arange(6, dtype=np.int64) + 100
            self.flux = np.arange(20*2*3*4, dtype=np.float32).reshape((20, 2, 3, 4))
            self.velocity = np.arange(20*3, dtype=np.float32).reshape((20, 3))
            self.labels = np.array(["a", "b", "c"])
            self.bad_time = np.arange(8)
            self.DQF = np.zeros(20, dtype=np.int16)

        def cdf_info(self):
            return SimpleNamespace(zVariables=[
                "Epoch", "Epoch_ALT", "flux", "velocity",
                "labels", "bad_time", "DQF"], rVariables=[])

        def varinq(self, name):
            arr = self.varget(name)
            return SimpleNamespace(Rec_Vary=name != "labels",
                                   Last_Rec=len(arr)-1, Dim_Sizes=list(arr.shape[1:]))

        def varattsget(self, name):
            dep = ("Epoch_ALT" if name == "bad_time_alt" else
                   "Epoch" if name not in ("labels", "bad_time") else "")
            return {"DEPEND_0": dep, "UNITS": "test", "FILLVAL": -1e31}

        def varget(self, name, *, startrec=None, endrec=None):
            arr = getattr(self, name)
            return arr if startrec is None else arr[startrec:endrec+1]

    fake = ModuleType("cdflib")
    fake.cdfepoch = SimpleNamespace(to_datetime=lambda values:
        np.datetime64("2018-11-06T00:00:00", "ns") +
        np.asarray(values).astype("timedelta64[s]"))
    with tempfile.TemporaryDirectory() as directory, patch.dict(sys.modules, {"cdflib": fake}):
        csv_path = Path(directory) / "fields.csv"
        pd.DataFrame([{"event_id": "../../not_a_filepath",
                       "event_start_utc": "2018-11-06T00:00:03Z",
                       "event_end_utc": "2018-11-06T00:00:06Z",
                       "event_end_exclusive_utc": "2018-11-06T00:00:07Z"}]).to_csv(
                           csv_path, index=False)
        event = load_event_intervals(csv_path, before_seconds=1,
                                     after_seconds=1)[0]
        assert event["window_start_utc"] == "2018-11-06T00:00:02+00:00"
        assert event["window_end_exclusive_utc"] == "2018-11-06T00:00:08+00:00"
        assert "/" not in event["event_key"]
        report, arrays = native_event_slice(FakeCDF(), event["window_start_utc"],
                                             event["window_end_exclusive_utc"])
        assert arrays["flux"].shape == (6, 2, 3, 4)
        assert np.array_equal(arrays["flux"], FakeCDF().flux[2:8])
        assert np.array_equal(arrays["velocity"], FakeCDF().velocity[2:8])
        assert arrays["Epoch"].tolist() == list(range(2, 8))
        assert arrays["labels"].tolist() == ["a", "b", "c"]
        assert "bad_time" not in arrays
        assert report["variables"]["bad_time"]["status"] == "skipped_unresolved_time_axis"
        assert report["time_axes"]["Epoch"]["records_in_window"] == 6
        assert report["time_axes"]["Epoch_ALT"]["records_in_window"] == 0
        archive = Path(directory) / "test.npz"
        np.savez_compressed(archive, **arrays)
        with np.load(archive, allow_pickle=False) as saved:
            assert saved["flux"].shape == (6, 2, 3, 4)
            assert np.array_equal(saved["Epoch"], np.arange(2, 8))
        segments = list(split_utc_days("2018-11-05T23:59:58Z",
                                       "2018-11-06T00:00:02Z"))
        assert len(segments) == 2 and segments[0][2] == segments[1][1]
        # Full offline pipeline: archive + manifest + resumability, no NASA I/O.
        fake.CDF = lambda path: FakeCDF()
        with patch.object(sys.modules[__name__], "event_source",
                          return_value={"status": "cached",
                                        "path": str(Path(directory) / "fixture.cdf"),
                                        "mode": "synthetic_fixture"}):
            target = Path(directory) / "event_out"
            first = extract_events(csv_path, target, Path(directory) / "cache",
                                   ["spc-l3i"], max_events=1,
                                   before_seconds=1, after_seconds=1)
            assert first["status_counts"] == {"saved": 1}
            saved_files = list(target.glob("spc-l3i/*.npz"))
            assert len(saved_files) == 1
            with np.load(saved_files[0], allow_pickle=False) as saved:
                assert saved["flux"].shape == (6, 2, 3, 4)
            assert pd.read_csv(target / "event_manifest.csv").iloc[0]["records_primary_epoch"] == 6
            again = extract_events(csv_path, target, Path(directory) / "cache",
                                   ["spc-l3i"], max_events=1,
                                   before_seconds=1, after_seconds=1)
            assert again["status_counts"] == {"resumed": 1}
    print("PASS: FIELDS exact UTC event bounds, flanks, 3D arrays, flags, "
          "time axes, midnight split, full NPZ/manifest pipeline and resume")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--list-products", action="store_true")
    p.add_argument("--extract-events", action="store_true",
                   help="Extract all available native SWEAP variables during exact FIELDS event intervals")
    p.add_argument("--events-csv", type=Path, help="FIELDS CSV with event_start_utc and exclusive UTC end")
    p.add_argument("--event-products", default="all",
                   help="Comma-separated SWEAP keys; default all 13 products")
    p.add_argument("--event-out", type=Path, default=Path("data/sweap/events"),
                   help="NPZ archives, per-segment JSON and coverage manifest")
    p.add_argument("--max-events", type=int, default=0,
                   help="Only process the first N catalog events; 0 = all 599")
    p.add_argument("--before-seconds", type=float, default=10,
                   help="Additional seconds preceding the FIELDS event")
    p.add_argument("--after-seconds", type=float, default=10,
                   help="Additional seconds following the FIELDS event")
    p.add_argument("--max-variable-mib", type=float, default=128,
                   help="Maximum uncompressed bytes per individual event variable")
    p.add_argument("--no-resume", dest="resume", action="store_false",
                   help="Re-extract existing event NPZ archives")
    p.add_argument("--survey", action="store_true",
                   help="Inspect all 8 SWEAP measurement families / 13 public NASA product variants")
    p.add_argument("--catalog-only", action="store_true",
                   help="Survey NASA file availability without downloading any CDF")
    p.add_argument("--max-file-mib", type=float, default=50,
                   help="Survey download cap per CDF file in MiB (default 50)")
    p.add_argument("--max-total-mib", type=float, default=250,
                   help="Survey total download cap in MiB (default 250)")
    p.add_argument("--subset-on-limit", action="store_true",
                   help="Ask NASA for an 8-minute CDF subset when the daily file exceeds the size cap")
    p.add_argument("--families", help="Optional comma-separated survey families (see --list-products)")

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
        test()
        self_test()
        event_self_test()
        return
    if args.list_products:
        for key, (family, dataset, _) in CATALOG.items():
            print(f"{key:25} {family:22} {dataset}")
        return
    if args.extract_events:
        if args.survey:
            p.error("Choose --survey or --extract-events, not both")
        if args.events_csv is None or not args.events_csv.is_file():
            p.error("--extract-events requires an existing --events-csv FIELDS catalog")
        if args.max_events < 0:
            p.error("--max-events must be >= 0 (0 means all events)")
        if any(not np.isfinite(v) for v in
               (args.before_seconds, args.after_seconds,
                args.max_file_mib, args.max_total_mib, args.max_variable_mib)):
            p.error("Event time and size limits must be finite numbers")
        if args.before_seconds < 0 or args.after_seconds < 0:
            p.error("Event flanks must be nonnegative")
        if min(args.max_file_mib, args.max_total_mib, args.max_variable_mib) <= 0:
            p.error("Event size limits must be positive")
        event_products = (list(CATALOG) if args.event_products == "all"
                          else [x.strip() for x in args.event_products.split(",") if x.strip()])
        unknown = sorted(set(event_products) - set(CATALOG))
        if unknown or not event_products:
            p.error("Invalid --event-products: " + ", ".join(unknown or ["empty selection"]))
        # Do not duplicate repeated product keys within one run.
        event_products = list(dict.fromkeys(event_products))
        extract_events(args.events_csv, args.event_out, args.cache_dir, event_products,
                       max_events=args.max_events, before_seconds=args.before_seconds,
                       after_seconds=args.after_seconds,
                       max_file_mib=args.max_file_mib, max_total_mib=args.max_total_mib,
                       max_variable_mib=args.max_variable_mib, resume=args.resume)
        return
    if args.survey:
        if not args.date:
            p.error("--survey requires --date YYYY-MM-DD")
        if any(not np.isfinite(v) or v <= 0 for v in
               (args.minutes, args.max_file_mib, args.max_total_mib)):
            p.error("--minutes, --max-file-mib and --max-total-mib must be positive")
        all_families = {spec[0] for spec in CATALOG.values()}
        wanted = {s.strip() for s in args.families.split(",")} if args.families else None
        if wanted and (wanted - all_families):
            p.error("Unknown --families: " + ", ".join(sorted(wanted - all_families)))
        day = date.fromisoformat(args.date)
        output = args.json or args.cache_dir / f"survey_{day:%Y%m%d}.json"
        survey(day, args.cache_dir, output, minutes=args.minutes, start=args.start,
               catalog_only=args.catalog_only, max_file_mib=args.max_file_mib,
               max_total_mib=args.max_total_mib, families=wanted,
               subset_on_limit=args.subset_on_limit)
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
