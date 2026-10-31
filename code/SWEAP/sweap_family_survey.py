"""Preliminary survey of 8 SWEAP measurement families (13 NASA CDF variants).

This tool inspects real CDF metadata and a few native records, without
resampling or interpreting measurements as science-quality plasma moments.
Import through sweap_data_reader.py --survey.
"""
from datetime import date, timedelta
from pathlib import Path
from urllib.parse import urlparse
import json
import re

import numpy as np
import pandas as pd
import requests

CDAS = "https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/sp_phys/datasets"

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
