"""Parker Solar Probe FIELDS magnetic-field reader.

One interface for the *two distinct NASA data products*:
  get_day(date, resolution="1min")  -- one-minute Level-2 RTN averages
  get_day(date, resolution="4cyc")  -- four-samples-per-cycle Level-2 RTN

Both return a UTC-indexed DataFrame with BR, BT, BN and Bmag (nT).
The products are not interchangeable: the low-cadence product is not
interpolated into the high-cadence measurements. Event detection is handled
separately by reference_detector_v2.py.
"""

import argparse
from datetime import date, timedelta
from pathlib import Path
from urllib.parse import urlparse
import re

import cdflib
import numpy as np
import pandas as pd
import requests

BASE_API = "https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/sp_phys/datasets"
PRODUCTS = {
    "1min": {
        "dataset": "PSP_FLD_L2_MAG_RTN_1MIN",
        "vector": "psp_fld_l2_mag_RTN_1min",
        "epoch": "epoch_mag_RTN_1min",
        "filename_tag": "mag_rtn_1min",
        "cache": Path("data/fields"),
    },
    "4cyc": {
        "dataset": "PSP_FLD_L2_MAG_RTN_4_SA_PER_CYC",
        "vector": "psp_fld_l2_mag_RTN_4_Sa_per_Cyc",
        "epoch": "epoch_mag_RTN_4_Sa_per_Cyc",
        "filename_tag": "mag_rtn_4_sa_per_cyc",
        "cache": Path("data/fields_highres"),
    },
}
# Instrument quality flags: thruster firing, spacecraft roll and calibration.
BAD_BITS = 2 | 8 | 16


def days(start, n):
    first = date.fromisoformat(start)
    return [first + timedelta(days=i) for i in range(n)]


def version(name):
    match = re.search(r"_v(\d+)\.cdf$", name, re.I)
    return int(match.group(1)) if match else -1


def _product(resolution):
    if resolution not in PRODUCTS:
        raise ValueError(f"Unknown resolution {resolution!r}; choose 1min or 4cyc")
    return PRODUCTS[resolution]


def download_day(day, resolution="1min"):
    """Find a cached CDF or download the requested official NASA product."""
    product = _product(resolution)
    day = pd.Timestamp(day).date()
    marker = day.strftime("%Y%m%d")
    cache = product["cache"]
    # Preserve the search locations used by the pre-merge high-cadence reader.
    folders = (cache,) if resolution == "1min" else (
        cache, Path("data/fields"), Path("data")
    )
    cached = [path for folder in folders for path in folder.glob("*.cdf")
              if product["filename_tag"] in path.name.lower()
              and marker in path.name]
    if cached:
        return max(cached, key=lambda path: version(path.name))

    cache.mkdir(parents=True, exist_ok=True)
    stop = day + timedelta(days=1)
    api = f'{BASE_API}/{product["dataset"]}/orig_data'
    url = f"{api}/{day:%Y%m%d}T000000Z,{stop:%Y%m%d}T000000Z/"
    response = requests.get(url, headers={"Accept": "application/json"}, timeout=90)
    response.raise_for_status()
    items = response.json().get("FileDescription", [])
    if isinstance(items, dict):
        items = [items]
    matches = [(Path(urlparse(item["Name"]).path).name, item["Name"])
               for item in items if item.get("Name", "").lower().endswith(".cdf")
               and marker in item["Name"]]
    if not matches:
        raise FileNotFoundError(f"No FIELDS {resolution} CDF found for {day}")
    filename, remote = max(matches, key=lambda item: version(item[0]))
    target = cache / filename
    partial = cache / (filename + ".part")
    try:
        with requests.get(remote, stream=True, timeout=(30, 180)) as stream:
            stream.raise_for_status()
            with partial.open("wb") as out:
                for chunk in stream.iter_content(1024 * 1024):
                    if chunk:
                        out.write(chunk)
        if partial.stat().st_size < 1024:
            raise ValueError("Downloaded CDF is unexpectedly small")
        partial.replace(target)
    finally:
        partial.unlink(missing_ok=True)
    return target


def get_day(day, resolution="1min"):
    """Return calibrated RTN components and |B| for one day in UTC.

    The 4cyc product is the input for the detector's later 1-second resampling;
    the 1min product provides its long-time background estimate.
    """
    product = _product(resolution)
    cdf = cdflib.CDF(str(download_day(day, resolution=resolution)))
    names = {name.lower(): name for name in (
        list(cdf.cdf_info().zVariables) + list(cdf.cdf_info().rVariables)
    )}
    key = product["vector"].lower()
    epoch_key = product["epoch"].lower()
    if key not in names or epoch_key not in names:
        raise ValueError(f"FIELDS {resolution}: expected RTN magnetic vector/epoch missing")
    vector_name = names[key]
    magnetic = np.asarray(cdf.varget(vector_name), dtype=float)
    timestamps = pd.to_datetime(cdflib.cdfepoch.to_datetime(cdf.varget(names[epoch_key])))
    if magnetic.ndim != 2 or magnetic.shape[1] != 3 or len(magnetic) != len(timestamps):
        raise ValueError(f"FIELDS {resolution}: expected an N x 3 RTN array")
    fill = cdf.varattsget(vector_name).get("FILLVAL")
    if fill is not None and np.isscalar(fill):
        magnetic[magnetic == float(fill)] = np.nan
    magnetic[(~np.isfinite(magnetic)) | (np.abs(magnetic) > 1e6)] = np.nan
    frame = pd.DataFrame(magnetic, index=timestamps, columns=["BR", "BT", "BN"])
    frame = frame[~frame.index.duplicated(keep="last")].sort_index()

    if {"psp_fld_l2_quality_flags", "epoch_quality_flags"}.issubset(names):
        flag_dates = pd.to_datetime(cdflib.cdfepoch.to_datetime(
            cdf.varget(names["epoch_quality_flags"])
        ))
        flag_values = np.asarray(cdf.varget(names["psp_fld_l2_quality_flags"])).ravel()
        if len(flag_dates) != len(flag_values):
            raise ValueError(f"FIELDS {resolution}: invalid quality flag lengths")
        flags = pd.Series(flag_values, index=flag_dates).sort_index()
        flags = flags[~flags.index.duplicated(keep="last")]
        aligned = flags.reindex(frame.index, method="nearest",
                                tolerance=pd.Timedelta(seconds=35))
        bad = ((aligned.fillna(0).astype("int64").to_numpy() & BAD_BITS) != 0)
        # Preserve prior product-specific treatment of missing flag timestamps.
        # The 4cyc reader invalidated them; the original 1min reader did not.
        if resolution == "4cyc":
            bad |= aligned.isna().to_numpy()
        frame.loc[bad, ["BR", "BT", "BN"]] = np.nan

    frame["Bmag"] = np.sqrt((frame[["BR", "BT", "BN"]] ** 2).sum(
        axis=1, min_count=3
    ))
    frame.index.name = "UTC"
    return frame


def read_high(day):
    """Small compatibility wrapper for the detector; still one reader module."""
    return get_day(day, resolution="4cyc")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Read PSP FIELDS RTN magnetic data")
    parser.add_argument("--start", default="2018-11-06")
    parser.add_argument("--days", type=int, default=1)
    parser.add_argument("--resolution", choices=tuple(PRODUCTS), default="1min")
    args = parser.parse_args()
    for day in days(args.start, args.days):
        df = get_day(day, resolution=args.resolution)
        print(f"{day} [{args.resolution}]: {df.Bmag.notna().sum()} "
              f"valid samples of {len(df)}")
