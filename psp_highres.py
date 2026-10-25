"""Load calibrated PSP/FIELDS L2 RTN magnetic-field vectors at 4 samples/cycle.

Uses CDAWeb's PSP_FLD_L2_MAG_RTN_4_SA_PER_CYC product and caches CDF files
locally in data/fields_highres/.  This is a newly consolidated loader;
confirm equivalence to older local survey routines before reproducing
previously published internal catalogs.
"""
from datetime import timedelta
from functools import lru_cache
from pathlib import Path
from urllib.parse import urlparse
import re

import cdflib
import numpy as np
import pandas as pd
import requests

CACHE = Path("data/fields_highres")
BASE = ("https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/sp_phys/"
        "datasets/PSP_FLD_L2_MAG_RTN_4_SA_PER_CYC/orig_data")
VECTOR = "psp_fld_l2_mag_RTN_4_Sa_per_Cyc"
EPOCH = "epoch_mag_RTN_4_Sa_per_Cyc"
BAD_BITS = 2 | 8 | 16  # Thruster, spacecraft roll, MAG calibration


def _version(path):
    m = re.search(r"_v(\d+)\.cdf$", path.name, re.I)
    return int(m.group(1)) if m else -1


def _file(day):
    marker = day.strftime("%Y%m%d")
    cached = [p for folder in (CACHE, Path("data/fields"), Path("data"))
              for p in folder.glob("*.cdf")
              if "mag_rtn_4_sa_per_cyc" in p.name.lower() and marker in p.name]
    if cached:
        return max(cached, key=_version)
    CACHE.mkdir(parents=True, exist_ok=True)
    stop = day + timedelta(days=1)
    url = f"{BASE}/{day:%Y%m%d}T000000Z,{stop:%Y%m%d}T000000Z/"
    response = requests.get(url, headers={"Accept": "application/json"}, timeout=90)
    response.raise_for_status()
    descriptions = response.json().get("FileDescription", [])
    if isinstance(descriptions, dict):
        descriptions = [descriptions]
    found = [(Path(urlparse(f["Name"]).path).name, f["Name"])
             for f in descriptions if f.get("Name", "").lower().endswith(".cdf")
             and marker in f["Name"]]
    if not found:
        raise FileNotFoundError(f"No 4cyc FIELDS CDF available for {day}")
    name, remote = max(found, key=lambda item: _version(Path(item[0])))
    target = CACHE / name
    temporary = target.with_suffix(target.suffix + ".part")
    try:
        with requests.get(remote, stream=True, timeout=(30, 180)) as stream:
            stream.raise_for_status()
            with temporary.open("wb") as out:
                for chunk in stream.iter_content(1024 * 1024):
                    if chunk:
                        out.write(chunk)
        if temporary.stat().st_size < 1024:
            raise ValueError("Downloaded CDF is too small")
        temporary.replace(target)
    finally:
        temporary.unlink(missing_ok=True)
    return target


@lru_cache(maxsize=3)
def read_high(day):
    """Return BR, BT, BN and |B| in nT, indexed by UTC timestamps."""
    day = pd.Timestamp(day).date()
    cdf = cdflib.CDF(str(_file(day)))
    names = {n.lower(): n for n in
             list(cdf.cdf_info().zVariables) + list(cdf.cdf_info().rVariables)}
    if VECTOR.lower() not in names or EPOCH.lower() not in names:
        raise ValueError(f"Missing expected RTN vector or epoch: {list(names)}")
    var = names[VECTOR.lower()]
    magnetic = np.asarray(cdf.varget(var), dtype=float)
    epochs = pd.to_datetime(cdflib.cdfepoch.to_datetime(cdf.varget(names[EPOCH.lower()])))
    if magnetic.ndim != 2 or magnetic.shape[1] != 3 or len(magnetic) != len(epochs):
        raise ValueError("Unexpected RTN vector shape")
    fill = cdf.varattsget(var).get("FILLVAL")
    if fill is not None and np.isscalar(fill):
        magnetic[magnetic == float(fill)] = np.nan
    magnetic[(~np.isfinite(magnetic)) | (np.abs(magnetic) > 1e6)] = np.nan
    df = pd.DataFrame(magnetic, index=epochs, columns=["BR", "BT", "BN"])
    df = df[~df.index.duplicated(keep="last")].sort_index()
    if {"psp_fld_l2_quality_flags", "epoch_quality_flags"}.issubset(names):
        flag_times = pd.to_datetime(cdflib.cdfepoch.to_datetime(
            cdf.varget(names["epoch_quality_flags"])))
        flags = pd.Series(np.asarray(cdf.varget(names["psp_fld_l2_quality_flags"])).ravel(),
                          index=flag_times).sort_index()
        flags = flags[~flags.index.duplicated(keep="last")]
        aligned = flags.reindex(df.index, method="nearest",
                                tolerance=pd.Timedelta(seconds=35))
        invalid = aligned.isna().to_numpy() | (
            (aligned.fillna(0).astype("int64").to_numpy() & BAD_BITS) != 0)
        df.loc[invalid, ["BR", "BT", "BN"]] = np.nan
    df["Bmag"] = np.sqrt((df[["BR", "BT", "BN"]] ** 2).sum(axis=1, min_count=3))
    df.index.name = "UTC"
    return df
