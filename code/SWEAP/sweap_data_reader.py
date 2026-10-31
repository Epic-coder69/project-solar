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


EVENT_EXTRACTOR_VERSION = "native_sweap_v2_context_quality_audit"

def cdf_json(value):
    """Represent standard CDF metadata scalars/arrays in JSON."""
    if isinstance(value, dict):
        return {str(k): cdf_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [cdf_json(v) for v in value]
    if isinstance(value, np.ndarray):
        return cdf_json(value.tolist())
    if isinstance(value, np.generic):
        return cdf_json(value.item())
    if isinstance(value, (pd.Timestamp, date)):
        return value.isoformat()
    if isinstance(value, (bytes, bytearray)):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, float) and not np.isfinite(value):
        return str(value)
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)

def count_valid_native(values, attrs):
    """Non-destructive validity counts, not a replacement for quality flags."""
    a = np.asarray(values)
    if not np.issubdtype(a.dtype, np.number):
        return None
    good = np.isfinite(a)
    for name, comparison in (
        ("FILLVAL", "ne"), ("VALIDMIN", "ge"), ("VALIDMAX", "le")):
        if name not in attrs:
            continue
        try:
            bound = np.asarray(attrs[name], dtype=float)
            if comparison == "ne":
                good &= a != bound
            elif comparison == "ge":
                good &= a >= bound
            else:
                good &= a <= bound
        except (TypeError, ValueError):
            pass
    return {"total_elements": int(a.size),
            "not_fill_and_in_valid_range": int(good.sum()),
            "nonfinite_or_fill_or_out_of_range": int((~good).sum()),
            "rows_with_any_valid_element": int(
                good.reshape((len(a), -1)).any(axis=1).sum())
            if a.ndim and len(a) else None}

def epoch_phase_counts(times, event_start, event_end, before_seconds, after_seconds):
    """Classify EPOCH tags, not integration intervals."""
    start, stop = utc(event_start), utc(event_end)
    pre = start - pd.Timedelta(seconds=before_seconds)
    post = stop + pd.Timedelta(seconds=after_seconds)
    return {
        "during_event": int(np.sum((times >= start) & (times < stop))),
        "immediate_pre": int(np.sum((times >= pre) & (times < start))),
        "immediate_post": int(np.sum((times >= stop) & (times < post))),
        "outer_support_before": int(np.sum(times < pre)),
        "outer_support_after": int(np.sum(times >= post))}


# Scientific observables promised by published product descriptions, but
# sometimes absent from an individual NASA-generated subset or source file.
# An absence is an explicit warning, not evidence that the native source has
# the same omission or permission to synthesize the missing measurements.
RECOMMENDED_SOURCE_VARIABLES = {
    "spc-l2i": ("measurement_time", "mv_hi"),
    "spa-electron-3d": ("THETA", "PHI", "QUALITY_FLAG"),
    "spi-alpha-l2": ("THETA", "PHI"),
}


def audit_native_quality(item):
    """Inspect saved native flags at FIELDS event EPOCHs, without changing data.

    A zero integer flag is provisionally 'unflagged', not a guarantee that
    viewing geometry, contamination or statistical quality is acceptable.
    """
    import cdflib

    archive = Path(item.get("archive", ""))
    if not archive.is_file():
        return {"status": "archive_not_found", "reason": str(archive)}
    try:
        with np.load(archive, allow_pickle=False) as saved:
            names = set(saved.files)
            epoch_name = next((v for v in names if v.lower() == "epoch"), None)
            if epoch_name is None:
                return {"status": "missing_epoch"}
            times = pd.DatetimeIndex(pd.to_datetime(
                cdflib.cdfepoch.to_datetime(saved[epoch_name]), utc=True))
            during = ((times >= utc(item["event_start_utc"])) &
                      (times < utc(item["event_end_exclusive_utc"])))
            flag_name = next((n for n in ("general_flag", "QUALITY_FLAG", "quality_flag", "DQF")
                              if n in names), None)
            if flag_name is None:
                return {"status": "no_flag_in_product", "event_epoch_count": int(during.sum())}
            flags = np.asarray(saved[flag_name])
            if len(flags) != len(times):
                return {"status": "quality_epoch_length_mismatch",
                        "flag": flag_name, "records": int(len(flags)),
                        "epoch_records": int(len(times))}
            values = flags[during]
            if not np.issubdtype(values.dtype, np.number):
                return {"status": "quality_non_numeric", "flag": flag_name}
            if values.ndim == 1:
                finite = np.isfinite(values)
                clean = finite & (values == 0)
                bad = finite & (values != 0)
                return {
                    "status": "audited", "flag": flag_name,
                    "event_epoch_count": int(len(values)),
                    "unflagged_zero": int(clean.sum()),
                    "flagged_nonzero": int(bad.sum()),
                    "nonfinite": int((~finite).sum()),
                    "unflagged_fraction": float(clean.mean()) if len(clean) else None,
                    "distinct_flag_values": {
                        str(v): int(n) for v, n in zip(*np.unique(values, return_counts=True))
                    }}
            # SPC DQF array [time, 32] is multidimensional; list nonzero
            # counts per bit separately, rather than collapsing to a badness.
            valid = np.isfinite(values)
            nonzero = valid & (values != 0)
            return {
                "status": "audited_multidimensional", "flag": flag_name,
                "event_epoch_count": int(len(values)),
                "shape": list(values.shape),
                "rows_with_any_nonzero_flag": int(nonzero.reshape((len(values), -1)).any(axis=1).sum()),
                "flag_nonzero_count_per_component": nonzero.sum(axis=0).tolist(),
                "nonfinite_elements": int((~valid).sum())}
    except (OSError, ValueError, TypeError, KeyError, IndexError) as exc:
        return {"status": "quality_read_error", "reason": str(exc)}


def audit_event_archives(directory):
    """Summarize genuine event-time EPOCHs and variables actually retained."""
    metadata = sorted(directory.glob("*/*.json"))
    rows = []
    for path in metadata:
        try:
            item = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        timing = item.get("variable_details", {}).get("time_axes", {})
        primary = next((x for x in timing if x.lower() == "epoch"), None)
        counts = timing.get(primary, {}).get("phase_counts", {}) if primary else {}
        quality = audit_native_quality(item)
        named = set(item.get("variable_details", {}).get("variables", {}))
        missing = [v for v in RECOMMENDED_SOURCE_VARIABLES.get(item.get("product"), ())
                   if v.lower() not in {n.lower() for n in named}]
        rows.append({
            "event_id": item.get("event_id"), "product": item.get("product"),
            "status": item.get("status"),
            "during_epoch": counts.get("during_event"),
            "immediate_pre_epoch": counts.get("immediate_pre"),
            "immediate_post_epoch": counts.get("immediate_post"),
            "support_before_epoch": counts.get("outer_support_before"),
            "support_after_epoch": counts.get("outer_support_after"),
            "saved_variables": item.get("saved_variables"),
            "not_saved_variables": item.get("not_saved_variables"),
            "read_errors": item.get("read_errors"),
            "science_variables_with_event_values": item.get("science_variables_with_event_values"),
            "all_numeric_variables_with_event_values": item.get("all_numeric_variables_with_event_values"),
            "spc_steps_during_event": item.get("spc_steps_during_event"),
            "quality_flag_name": quality.get("flag"),
            "quality_audit_status": quality.get("status"),
            "event_unflagged_zero": quality.get("unflagged_zero"),
            "event_flagged_nonzero": quality.get("flagged_nonzero"),
            "event_quality_flag_fraction": quality.get("unflagged_fraction"),
            "missing_recommended_source_variables": ";".join(missing),
            "recommended_source_variable_gap_count": len(missing),
            "source": item.get("source_file"),
            "extractor_version": item.get("extractor_version"),
            "native_flag_audit": quality})
    if not rows:
        raise FileNotFoundError(f"No SWEAP product metadata JSON found under {directory}")
    table = pd.DataFrame(rows)
    # Native flag component counts and distinct values are serialized in JSON,
    # while the CSV remains simple enough to read in any spreadsheet.
    table.assign(native_flag_audit=table.native_flag_audit.apply(
        lambda q: json.dumps(q, default=cdf_json))).to_csv(
            directory / "coverage_audit.csv", index=False)
    report = {
        "product_segments": len(table),
        "events": int(table.event_id.nunique()),
        "status_counts": table.status.value_counts().to_dict(),
        "segments_with_during_epoch": int((table.during_epoch.fillna(0) > 0).sum()),
        "segments_with_zero_during_epoch": int((table.during_epoch.fillna(0) == 0).sum()),
        "segments_with_skipped_variables": int((table.not_saved_variables.fillna(0) > 0).sum()),
        "segments_with_read_errors": int((table.read_errors.fillna(0) > 0).sum()),
        "segments_with_valid_science_values_during_event": int(
            (table.science_variables_with_event_values.fillna(0) > 0).sum()),
        "segments_with_any_nonfill_numeric_values_during_event": int(
            (table.all_numeric_variables_with_event_values.fillna(0) > 0).sum()),
        "segments_with_missing_recommended_source_variables": int(
            (table.recommended_source_variable_gap_count > 0).sum()),
        "missing_source_variables_by_product": {
            str(row["product"]): row["missing_recommended_source_variables"].split(";")
            for _, row in table.iterrows() if row["recommended_source_variable_gap_count"] > 0
        },
        "event_native_quality_by_event_product": {
            f'{row["event_id"]}|{row["product"]}|{row["source"]}':
                row["native_flag_audit"]
            for _, row in table.iterrows()
        },
        "segments_with_unflagged_zero": int(
            (table.event_unflagged_zero.fillna(0) > 0).sum()),
        "notes": [
            "EPOCH tag counts are not exposure-overlap counts.",
            "A nonzero EPOCH count is not a science-quality approval.",
            "Zero quality flags are provisional only; apply instrument-specific flag/FOV definitions.",
            "Missing recommended source variables may result from NASA subsets or product generation.",
            "Refer to per-event product JSON for exact variables and NASA caveats."]}
    (directory / "coverage_audit_summary.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("[SWEAP COVERAGE AUDIT]", report, flush=True)
    return report


def bundle_event_audit(directory):
    """Share small event metadata and coverage audit, never huge native NPZs."""
    from zipfile import ZipFile, ZIP_DEFLATED
    archive = directory / "sweap_audit_bundle.zip"
    samples = sorted(directory.glob("*/*.json"))
    for name in ("event_manifest.csv", "summary.json",
                 "coverage_audit.csv", "coverage_audit_summary.json"):
        candidate = directory / name
        if candidate.is_file():
            samples.append(candidate)
    if not samples:
        raise FileNotFoundError(f"No event metadata to bundle under {directory}")
    with ZipFile(archive, "w", compression=ZIP_DEFLATED, compresslevel=6) as zipped:
        for source in samples:
            zipped.write(source, arcname=source.relative_to(directory))
    print(f"[SWEAP AUDIT BUNDLE] {archive} ({len(samples)} metadata files)", flush=True)
    return archive


def load_event_intervals(path, max_events=0, before_seconds=10, after_seconds=10, support_seconds=120):
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
        begin = start - pd.Timedelta(seconds=before_seconds + support_seconds)
        stop = end + pd.Timedelta(seconds=after_seconds + support_seconds)
        # No user-provided event ID appears in output paths (path traversal safe).
        key = hashlib.sha256(f"{event_id}|{start.isoformat()}|{end.isoformat()}".encode()).hexdigest()[:16]
        result.append({"event_id": event_id, "event_key": key,
                       "event_start_utc": start.isoformat(),
                       "event_end_exclusive_utc": end.isoformat(),
                       "window_start_utc": begin.isoformat(),
                       "window_end_exclusive_utc": stop.isoformat(),
                       "before_seconds": float(before_seconds),
                       "after_seconds": float(after_seconds),
                       "support_seconds": float(support_seconds)})
    return result


def split_utc_days(begin, stop):
    """Split [begin,stop) at UTC midnight with no duplicated boundaries."""
    a, b = utc(begin), utc(stop)
    while a < b:
        next_stop = min(b, a.floor("D") + pd.Timedelta(days=1))
        yield a.date(), a, next_stop
        a = next_stop


def native_event_slice(cdf, begin, stop, max_variable_mib=0,
                       event_start=None, event_end=None,
                       before_seconds=10, after_seconds=10, product=None):
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
              "full_native_arrays": True, "extractor_version": EVENT_EXTRACTOR_VERSION,
              "metadata": {"global_attributes": cdf_json(cdf.globalattsget())
                           if hasattr(cdf, "globalattsget") else {}}}
    selections = {}
    left, right = utc(begin), utc(stop)
    event_start = utc(event_start) if event_start is not None else left
    event_end = utc(event_end) if event_end is not None else right
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
            "last_utc": times[ids[-1]].isoformat() if len(ids) else None,
            "phase_counts": epoch_phase_counts(times[ids], event_start,
                                                event_end, before_seconds,
                                                after_seconds)}
    arrays = {}
    max_bytes = int(max_variable_mib * 1048576) if max_variable_mib else None
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
            if max_bytes is not None and slab.nbytes > max_bytes:
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
            during_counts = None
            if varying and np.issubdtype(slab.dtype, np.number):
                source_times = timestamps[dep][rowids]
                during = (source_times >= event_start) & (source_times < event_end)
                during_counts = count_valid_native(slab[during], attrs)
            result["variables"][name] = {
                "status": "saved", "rec_vary": varying,
                "depend_0": dep if varying else None,
                "time_axis_inferred": inferred,
                "shape": list(slab.shape), "dtype": str(slab.dtype),
                "units": str(attrs.get("UNITS", "")),
                "catdesc": str(attrs.get("CATDESC", ""))[:200],
                "fillval": str(attrs.get("FILLVAL", "")),
                "bytes": int(slab.nbytes),
                "validity_summary": count_valid_native(slab, attrs),
                "validity_during_event": during_counts,
                "cdf_variable_type": str(attrs.get("VAR_TYPE", "unknown")),
                "attributes": cdf_json(attrs),
                "cdf_data_type": str(getattr(var, "Data_Type_Description",
                                               getattr(var, "Data_Type", ""))),
                "cdf_record_dimensions": cdf_json(getattr(var, "Dim_Sizes", []) or [])}
        except (ValueError, OSError, TypeError, IndexError, AttributeError) as exc:
            result["variables"][name] = {
                "status": "read_error", "reason": f"{type(exc).__name__}: {exc}"}
            result["errors"][name] = f"{type(exc).__name__}: {exc}"
    # SPC L2 stores EPOCH at the start of the spectrum, not each bin's
    # measurement time. MEASUREMENT_TIME supplies each voltage step offset
    # in TT2000 nanoseconds. Convert both to a derived UTC nanosecond matrix.
    # The original CDF epoch/offset arrays remain unmodified in the NPZ.
    if product == "spc-l2i":
        epoch_name = next((n for n in arrays if n.lower() == "epoch"), None)
        offset_name = next((n for n in arrays if n.lower() == "measurement_time"), None)
        step_report = {"status": "unavailable", "reason": "Epoch or MEASUREMENT_TIME absent"}
        if epoch_name and offset_name and epoch_name in selections:
            offsets = np.asarray(arrays[offset_name])
            utc_times = timestamps[epoch_name][selections[epoch_name]]
            if (np.issubdtype(offsets.dtype, np.number) and offsets.ndim == 2
                    and offsets.shape[0] == len(utc_times)):
                offsets_float = offsets.astype(np.float64)
                # Unphysical delays and CDF fills are excluded, not converted
                # to real event times. Keep raw offsets for independent review.
                valid_offsets = np.isfinite(offsets_float) & (offsets_float >= 0) & (
                    offsets_float < 300_000_000_000)
                clean_offsets = np.where(valid_offsets, offsets_float, 0)
                measured = (utc_times.as_unit("ns").asi8[:, None] +
                            np.rint(clean_offsets).astype(np.int64))
                measured[~valid_offsets] = np.iinfo(np.int64).min
                arrays["derived_SPC_step_UTC_ns"] = measured
                measured_utc = pd.to_datetime(
                    measured[valid_offsets], unit="ns", utc=True)
                step_report = {
                    "status": "derived", "source_epoch": epoch_name,
                    "offset_variable": offset_name, "offset_units": "ns",
                    "derived_variable": "derived_SPC_step_UTC_ns",
                    "derived_units": "Unix UTC nanoseconds",
                    "invalid_step_timestamp_sentinel": int(np.iinfo(np.int64).min),
                    "valid_step_count": int(valid_offsets.sum()),
                    "steps_during_event": int(np.sum(
                        (measured_utc >= event_start) & (measured_utc < event_end))),
                    "steps_immediate_pre": int(np.sum(
                        (measured_utc >= event_start - pd.Timedelta(seconds=before_seconds))
                        & (measured_utc < event_start))),
                    "steps_immediate_post": int(np.sum(
                        (measured_utc >= event_end) & (
                            measured_utc < event_end + pd.Timedelta(seconds=after_seconds)))),
                    "note": "Times are step timestamps, NOT a claim about step integration duration."}
            else:
                step_report = {"status": "unavailable",
                               "reason": "Inconsistent or nonnumeric MEASUREMENT_TIME dimensions"}
        result["spc_l2_per_voltage_step_timing"] = step_report
    result["saved_variables"] = sum(v["status"] == "saved"
                                    for v in result["variables"].values())
    result["not_saved_variables"] = len(names) - result["saved_variables"]
    result["any_records_in_window"] = any(len(v) for v in selections.values())
    result["unlimited_variable_size"] = max_bytes is None
    result["saved_record_varying_variables"] = sum(
        v.get("status") == "saved" and v.get("rec_vary", False)
        for v in result["variables"].values())
    result["complete_variable_extraction"] = (
        result["not_saved_variables"] == 0 and not result["errors"])
    result["numeric_time_varying_variables_with_event_values"] = sum(
        v.get("status") == "saved" and v.get("rec_vary", False) and
        (v.get("validity_during_event") or {}).get("not_fill_and_in_valid_range", 0) > 0
        for v in result["variables"].values())
    result["numeric_science_variables_with_event_values"] = sum(
        v.get("status") == "saved" and
        str(v.get("cdf_variable_type", "")).lower() == "data" and
        (v.get("validity_during_event") or {}).get("not_fill_and_in_valid_range", 0) > 0
        for v in result["variables"].values())
    return result, arrays



def plan_subset_groups(events, gap_seconds=30, max_minutes=12):
    """Reuse a NASA time subset across nearby event windows on the same day.

    Conservative grouping avoids one remote request per event. Intervals
    always contain each event's COMPLETE support window; no downsampling.
    """
    daily = {}
    for event in events:
        for day, begin, end in split_utc_days(event["window_start_utc"],
                                              event["window_end_exclusive_utc"]):
            daily.setdefault(day, []).append((begin, end))
    result = {}
    gap = pd.Timedelta(seconds=gap_seconds)
    max_span = pd.Timedelta(minutes=max_minutes)
    for day, windows in daily.items():
        groups = []
        for begin, end in sorted(windows):
            if groups and begin <= groups[-1][1] + gap and (
                    max(end, groups[-1][1]) - groups[-1][0] <= max_span):
                groups[-1] = (groups[-1][0], max(groups[-1][1], end))
            else:
                groups.append((begin, end))
        result[day] = groups
    return result

def event_source(product, day, start, stop, cache, state, budget, file_mib,
                 groups=None, subset_state=None):
    """Reuse a daily CDF when small, else request a NASA time subset."""
    key = (product, str(day))
    if key not in state or state[key]["status"] == "error":
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
    request_start, request_stop = start, stop
    for lo, hi in (groups or {}).get(day, []):
        if lo <= start and hi >= stop:
            request_start, request_stop = lo, hi
            break
    subset_key = (product, str(day), request_start.isoformat(), request_stop.isoformat())
    if subset_state is not None and subset_key in subset_state:
        remembered = subset_state[subset_key]
        if remembered["status"] in ("cached", "downloaded"):
            return remembered
        if request_start != start or request_stop != stop:
            # A large merged NASA request may exceed its size cap, while
            # the much smaller individual event can still fit.
            return event_source(product, day, start, stop, cache, state,
                                budget, file_mib, groups={}, subset_state=subset_state)
    if budget[0] <= 0:
        return {"status": "skipped_total_size_budget"}
    try:
        minutes = float((request_stop - request_start) / pd.Timedelta(minutes=1))
        subset = subset_for_window(product, day, cache, request_start.isoformat(), minutes)
        if subset["status"] == "found":
            status, used = bounded_download(subset["url"], subset["path"],
                                            int(file_mib * 1048576), budget[0])
            budget[0] -= used
            subset["status"] = status
        if subset_state is not None:
            subset_state[subset_key] = subset
        if subset["status"] not in ("cached", "downloaded") and (
                request_start != start or request_stop != stop):
            return event_source(product, day, start, stop, cache, state,
                                budget, file_mib, groups={}, subset_state=subset_state)
        return subset
    except (requests.RequestException, OSError, ValueError) as exc:
        if request_start != start or request_stop != stop:
            return event_source(product, day, start, stop, cache, state,
                                budget, file_mib, groups={}, subset_state=subset_state)
        return {"status": "error", "reason": str(exc)}


def extract_events(csv_path, output_dir, cache, products, *, max_events=0,
                   before_seconds=10, after_seconds=10,
                   max_file_mib=50, max_total_mib=2048,
                   max_variable_mib=0, support_seconds=120, resume=True):
    """Write all available native measurements plus an explicit coverage manifest."""
    import cdflib
    events = load_event_intervals(csv_path, max_events, before_seconds, after_seconds, support_seconds)
    output_dir.mkdir(parents=True, exist_ok=True)
    budget = [int(max_total_mib * 1048576)]
    daily_state = {}
    subset_state = {}
    subset_groups = plan_subset_groups(events)
    rows = []
    print(f"[SWEAP EVENTS] {len(events)} FIELDS candidates; {len(products)} NASA products; "
          f"flanks {before_seconds}/{after_seconds}s, coarse context {support_seconds}s each side", flush=True)
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
                       "archive": str(path), "metadata": str(meta_path),
                       "extractor_version": EVENT_EXTRACTOR_VERSION,
                       "unlimited_variable_size": max_variable_mib == 0}
                if resume and path.is_file() and meta_path.is_file():
                    try:
                        prior = json.loads(meta_path.read_text(encoding="utf-8"))
                        if (prior.get("status") == "saved"
                                and prior.get("extractor_version") == EVENT_EXTRACTOR_VERSION
                                and prior.get("window_start_utc") == event["window_start_utc"]
                                and prior.get("window_end_exclusive_utc") == event["window_end_exclusive_utc"]
                                and prior.get("unlimited_variable_size") == (max_variable_mib == 0)
                                and prior.get("before_seconds") == float(before_seconds)
                                and prior.get("after_seconds") == float(after_seconds)
                                and prior.get("support_seconds") == float(support_seconds)
                                and prior.get("source_file")
                                and path.stat().st_size > 100):
                            row.update({"status": "resumed",
                                        "saved_variables": prior.get("saved_variables", 0),
                                        "records_primary_epoch": prior.get("records_primary_epoch", 0),
                                        "records_during_event": prior.get("records_during_event", 0),
                                        "not_saved_variables": prior.get("not_saved_variables", 0),
                                        "complete_variable_extraction": prior.get("complete_variable_extraction", False),
                                        "science_variables_with_event_values": prior.get("science_variables_with_event_values", 0),
                                        "all_numeric_variables_with_event_values": prior.get("all_numeric_variables_with_event_values", 0),
                                        "read_errors": prior.get("read_errors", 0)})
                            rows.append(row)
                            continue
                    except (OSError, ValueError):
                        pass
                src = event_source(product, day, begin, end, cache, daily_state,
                                   budget, max_file_mib,
                                   groups=subset_groups, subset_state=subset_state)
                row["source_status"] = src["status"]
                row["source_file"] = src.get("path")
                row["source_url"] = src.get("url")
                row["source_kind"] = src.get("mode", "original_daily_cdf")
                if src["status"] not in ("cached", "downloaded"):
                    row.update(status="no_source", reason=src.get("reason", src["status"]))
                else:
                    try:
                        report, arrays = native_event_slice(
                            cdflib.CDF(str(src["path"])), begin, end, max_variable_mib,
                            event_start=event["event_start_utc"],
                            event_end=event["event_end_exclusive_utc"],
                            before_seconds=before_seconds, after_seconds=after_seconds,
                            product=product)
                        row["saved_variables"] = report["saved_variables"]
                        primary = next((n for n in report["time_axes"] if n.lower() == "epoch"),
                                       next(iter(report["time_axes"])))
                        row["records_primary_epoch"] = report["time_axes"][primary]["records_in_window"]
                        phases = report["time_axes"][primary]["phase_counts"]
                        row["records_during_event"] = phases["during_event"]
                        row["records_immediate_pre"] = phases["immediate_pre"]
                        row["records_immediate_post"] = phases["immediate_post"]
                        row["records_outer_support"] = (phases["outer_support_before"]
                                                        + phases["outer_support_after"])
                        row["complete_variable_extraction"] = report["complete_variable_extraction"]
                        row["science_variables_with_event_values"] = report[
                            "numeric_science_variables_with_event_values"]
                        row["all_numeric_variables_with_event_values"] = report[
                            "numeric_time_varying_variables_with_event_values"]
                        if product == "spc-l2i":
                            row["spc_steps_during_event"] = report.get(
                                "spc_l2_per_voltage_step_timing", {}).get("steps_during_event")
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
        "extractor_version": EVENT_EXTRACTOR_VERSION,
        "support_seconds": support_seconds, "max_variable_mib": max_variable_mib,
        "planned_grouped_time_subsets": sum(map(len, subset_groups.values())),
        "reused_grouped_subsets": len(subset_state),
        "product_segments_with_during_epoch": sum(r.get("records_during_event", 0) > 0 for r in rows),
        "product_segments_with_complete_variable_extraction": sum(
            r.get("complete_variable_extraction", False) for r in rows),
        "segments_with_valid_data_type_values_during_event": sum(
            r.get("science_variables_with_event_values", 0) > 0 for r in rows),
        "segments_with_nonfill_numeric_values_during_event": sum(
            r.get("all_numeric_variables_with_event_values", 0) > 0 for r in rows),
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
            self.MEASUREMENT_TIME = np.broadcast_to(
                np.array([0, 500_000_000, 900_000_000], dtype=np.int64),
                (20, 3)).copy()
            self.flux = np.arange(20*2*3*4, dtype=np.float32).reshape((20, 2, 3, 4))
            self.velocity = np.arange(20*3, dtype=np.float32).reshape((20, 3))
            self.labels = np.array(["a", "b", "c"])
            self.bad_time = np.arange(8)
            self.DQF = np.zeros(20, dtype=np.int16)

        def cdf_info(self):
            return SimpleNamespace(zVariables=[
                "Epoch", "Epoch_ALT", "MEASUREMENT_TIME", "flux", "velocity",
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
                                     after_seconds=1, support_seconds=0)[0]
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
        assert report["time_axes"]["Epoch"]["phase_counts"]["during_event"] == 6
        spc_report, spc_arrays = native_event_slice(
            FakeCDF(), event["window_start_utc"], event["window_end_exclusive_utc"],
            event_start=event["event_start_utc"],
            event_end=event["event_end_exclusive_utc"],
            before_seconds=1, after_seconds=1, product="spc-l2i")
        assert spc_arrays["derived_SPC_step_UTC_ns"].shape == (6, 3)
        assert spc_report["spc_l2_per_voltage_step_timing"]["steps_during_event"] == 12
        padded = load_event_intervals(csv_path, before_seconds=1,
                                      after_seconds=1, support_seconds=120)[0]
        assert utc(padded["window_start_utc"]) == utc(event["window_start_utc"]) - pd.Timedelta(seconds=120)
        grouped = plan_subset_groups([event, padded], gap_seconds=30)
        assert len(grouped[date(2018, 11, 6)]) == 1
        assert count_valid_native(np.array([1., -1e31, np.nan]),
                                  {"FILLVAL": -1e31})["not_fill_and_in_valid_range"] == 1
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
                                   before_seconds=1, after_seconds=1, support_seconds=0)
            assert first["status_counts"] == {"saved": 1}
            saved_files = list(target.glob("spc-l3i/*.npz"))
            assert len(saved_files) == 1
            with np.load(saved_files[0], allow_pickle=False) as saved:
                assert saved["flux"].shape == (6, 2, 3, 4)
            assert pd.read_csv(target / "event_manifest.csv").iloc[0]["records_primary_epoch"] == 6
            again = extract_events(csv_path, target, Path(directory) / "cache",
                                   ["spc-l3i"], max_events=1,
                                   before_seconds=1, after_seconds=1, support_seconds=0)
            assert again["status_counts"] == {"resumed": 1}
            audited = audit_event_archives(target)
            assert audited["segments_with_during_epoch"] == 1
            assert audited["segments_with_skipped_variables"] == 1
            assert len(audited["event_native_quality_by_event_product"]) == 1
            assert next(iter(audited["event_native_quality_by_event_product"].values()))["status"] == "audited"
            assert next(iter(audited["event_native_quality_by_event_product"].values()))["unflagged_zero"] == 4
            assert audited["segments_with_unflagged_zero"] == 1
    print("PASS: SWEAP event UTC windows, 120s context, grouped subsets, "
          "3D native arrays, fill masks, SPC nanosecond step times, "
          "metadata audit, NPZ/manifest and resume")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--list-products", action="store_true")
    p.add_argument("--bundle-audit", action="store_true",
                   help="Zip all SWEAP event JSON and coverage reports for sharing (no large NPZ data)")
    p.add_argument("--audit-events", action="store_true",
                   help="Audit saved SWEAP event JSON for event EPOCH counts and variable omissions")
    p.add_argument("--extract-events", action="store_true",
                   help="Extract all available native SWEAP variables during exact FIELDS event intervals")
    p.add_argument("--events-csv", type=Path, help="FIELDS CSV with event_start_utc and exclusive UTC end")
    p.add_argument("--event-products", default="all",
                   help="Comma-separated SWEAP keys; default all 13 products")
    p.add_argument("--event-out", type=Path, default=Path("data/sweap/events_v2"),
                   help="NPZ archives, per-segment JSON and coverage manifest")
    p.add_argument("--max-events", type=int, default=0,
                   help="Only process the first N catalog events; 0 = all 599")
    p.add_argument("--before-seconds", type=float, default=10,
                   help="Additional seconds preceding the FIELDS event")
    p.add_argument("--after-seconds", type=float, default=10,
                   help="Additional seconds following the FIELDS event")
    p.add_argument("--support-seconds", type=float, default=120,
                   help="Extra coarse-measurement context on each side of event and flanks (default 120s)")
    p.add_argument("--max-variable-mib", type=float, default=0,
                   help="Maximum uncompressed MiB per variable; 0 means no size cap (default)")
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
    if args.audit_events:
        if args.extract_events or args.survey:
            p.error("--audit-events cannot be combined with extract/survey")
        audit_event_archives(args.event_out)
        if args.bundle_audit:
            bundle_event_audit(args.event_out)
        return
    if args.bundle_audit:
        p.error("--bundle-audit must be combined with --audit-events")
    if args.extract_events:
        if args.survey:
            p.error("Choose --survey or --extract-events, not both")
        if args.events_csv is None or not args.events_csv.is_file():
            p.error("--extract-events requires an existing --events-csv FIELDS catalog")
        if args.max_events < 0:
            p.error("--max-events must be >= 0 (0 means all events)")
        if any(not np.isfinite(v) for v in
               (args.before_seconds, args.after_seconds, args.support_seconds,
                args.max_file_mib, args.max_total_mib, args.max_variable_mib)):
            p.error("Event time and size limits must be finite numbers")
        if args.before_seconds < 0 or args.after_seconds < 0 or args.support_seconds < 0:
            p.error("Event flanks must be nonnegative")
        if min(args.max_file_mib, args.max_total_mib) <= 0 or args.max_variable_mib < 0:
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
                       max_variable_mib=args.max_variable_mib,
                       support_seconds=args.support_seconds, resume=args.resume)
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
