"""Project Solar: PSP FIELDS one-minute magnetic-field measurements."""

import argparse
from datetime import date, timedelta
from pathlib import Path
from urllib.parse import urlparse
import re

import cdflib
import numpy as np
import pandas as pd
import requests

CACHE = Path("data/fields")
API = (
    "https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/"
    "sp_phys/datasets/PSP_FLD_L2_MAG_RTN_1MIN/orig_data"
)
VECTOR = "psp_fld_l2_mag_RTN_1min"
EPOCH = "epoch_mag_RTN_1min"

# Thruster firings, magnetometer rolls, MAG calibration
BAD_BITS = 2 | 8 | 16


def days(start, n):
    first = date.fromisoformat(start)
    return [first + timedelta(days=i) for i in range(n)]


def version(name):
    m = re.search(r"_v(\d+)\.cdf$", name, re.I)
    return int(m.group(1)) if m else -1


def download_day(day):
    CACHE.mkdir(parents=True, exist_ok=True)
    marker = day.strftime("%Y%m%d")

    existing = [
        p for p in CACHE.glob("*.cdf")
        if "mag_rtn_1min" in p.name.lower()
        and marker in p.name
    ]
    if existing:
        print("Cached:", marker)
        return max(existing, key=lambda p: version(p.name))

    stop = day + timedelta(days=1)
    url = (
        f"{API}/{day:%Y%m%d}T000000Z,"
        f"{stop:%Y%m%d}T000000Z/"
    )

    response = requests.get(
        url, headers={"Accept": "application/json"}, timeout=90
    )
    response.raise_for_status()

    items = response.json().get("FileDescription", [])
    if isinstance(items, dict):
        items = [items]

    found = [
        (Path(urlparse(x["Name"]).path).name, x["Name"])
        for x in items
        if x.get("Name")
        and marker in x["Name"]
        and x["Name"].lower().endswith(".cdf")
    ]

    if not found:
        raise RuntimeError(f"No NASA CDF found for {day}")

    name, file_url = max(
        found, key=lambda item: version(item[0])
    )

    target = CACHE / name
    partial = CACHE / (name + ".part")

    try:
        with requests.get(
            file_url, stream=True, timeout=(30, 180)
        ) as response:
            response.raise_for_status()
            with partial.open("wb") as output:
                for chunk in response.iter_content(1024 * 1024):
                    if chunk:
                        output.write(chunk)

        if partial.stat().st_size < 1024:
            raise ValueError("Downloaded file is too small")

        partial.replace(target)
    finally:
        partial.unlink(missing_ok=True)

    print("Downloaded:", target)
    return target


def get_day(day):
    path = download_day(day)
    cdf = cdflib.CDF(str(path))

    field = np.asarray(cdf.varget(VECTOR), dtype=float)
    timestamps = pd.to_datetime(
        cdflib.cdfepoch.to_datetime(cdf.varget(EPOCH))
    )

    if field.ndim != 2 or field.shape[1] != 3:
        raise ValueError("Expected three magnetic components")

    fill = cdf.varattsget(VECTOR).get("FILLVAL")
    if fill is not None and np.isscalar(fill):
        field[field == float(fill)] = np.nan

    field[
        ~np.isfinite(field) | (np.abs(field) > 1e6)
    ] = np.nan

    df = pd.DataFrame(
        field, index=timestamps, columns=["BR", "BT", "BN"]
    )
    df = df[~df.index.duplicated()].sort_index()

    names = (
        cdf.cdf_info().zVariables
        + cdf.cdf_info().rVariables
    )

    if {
        "psp_fld_l2_quality_flags",
        "epoch_quality_flags"
    }.issubset(set(names)):

        quality_times = pd.to_datetime(
            cdflib.cdfepoch.to_datetime(
                cdf.varget("epoch_quality_flags")
            )
        )
        quality_values = np.asarray(
            cdf.varget("psp_fld_l2_quality_flags")
        ).ravel()

        flags = pd.Series(
            quality_values, index=quality_times
        ).sort_index()

        flags = flags[
            ~flags.index.duplicated(keep="last")
        ]

        nearest = flags.reindex(
            df.index,
            method="nearest",
            tolerance=pd.Timedelta(seconds=35)
        ).fillna(0)

        bad = (
            nearest.to_numpy(dtype=np.int64) & BAD_BITS
        ) != 0

        df.loc[bad, ["BR", "BT", "BN"]] = np.nan

    df["Bmag"] = np.sqrt(
        (df[["BR", "BT", "BN"]] ** 2).sum(
            axis=1, min_count=3
        )
    )

    df.index.name = "UTC"
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", default="2018-11-06")
    parser.add_argument("--days", type=int, default=1)
    args = parser.parse_args()

    for day in days(args.start, args.days):
        df = get_day(day)
        print(
            f"{day}: {df.Bmag.notna().sum()} "
            f"valid of {len(df)} samples"
        )
