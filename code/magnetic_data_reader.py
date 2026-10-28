"""Parker Solar Probe FIELDS magnetic-field reader.

One interface for the *two distinct NASA data products*:
  get_day(date, resolution="1min")  -- one-minute Level-2 RTN averages
  get_day(date, resolution="4cyc")  -- four-samples-per-cycle Level-2 RTN

Both return a UTC-indexed DataFrame with BR, BT, BN and Bmag (nT).
The products are not interchangeable: the low-cadence product is not
interpolated into the high-cadence measurements. Event detection is handled
separately by magnetic_data_analysis.py.
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



# ---------------------------------------------------------------------------
# Full-cadence reading and *preliminary* switchback screening.
# A full-cadence CDF is split into four 6-hour files per day. A 1-second
# screening representation is kept in RAM; the original full-cadence CDFs
# remain available for detailed event analysis.
# ---------------------------------------------------------------------------
FULL_DATASET = "PSP_FLD_L2_MAG_RTN"
FULL_VECTOR = "psp_fld_l2_mag_RTN"
FULL_EPOCH = "epoch_mag_RTN"
FULL_CACHE = Path("data/fields_full")
FULL_NAME = re.compile(r"^psp_fld_l2_mag_rtn_(\d{10})_v\d+\.cdf$", re.I)


def _full_day_files(day):
    """Return all full-cadence segments for a day, not just the latest file."""
    day = pd.Timestamp(day).date()
    marker = day.strftime("%Y%m%d")
    FULL_CACHE.mkdir(parents=True, exist_ok=True)

    # A valid complete day should have four 6-hour segments. Retain
    # downloaded segments even if CDAWeb temporarily becomes unavailable.
    known = {}
    for folder in (FULL_CACHE, Path("data/fields"), Path("data")):
        for path in folder.glob("*.cdf"):
            match = FULL_NAME.fullmatch(path.name)
            if match and match.group(1).startswith(marker):
                key = match.group(1)
                if key not in known or version(path.name) > version(known[key].name):
                    known[key] = path

    expected = {f"{marker}{hour:02d}" for hour in (0, 6, 12, 18)}
    if expected.issubset(known):
        return [known[key] for key in sorted(expected)]

    stop = day + timedelta(days=1)
    api = f"{BASE_API}/{FULL_DATASET}/orig_data"
    url = f"{api}/{day:%Y%m%d}T000000Z,{stop:%Y%m%d}T000000Z/"
    response = requests.get(url, headers={"Accept": "application/json"}, timeout=90)
    response.raise_for_status()
    descriptions = response.json().get("FileDescription", [])
    if isinstance(descriptions, dict):
        descriptions = [descriptions]
    remote = {}
    for entry in descriptions:
        address = entry.get("Name", "")
        filename = Path(urlparse(address).path).name
        match = FULL_NAME.fullmatch(filename)
        if not match or not match.group(1).startswith(marker):
            continue
        key = match.group(1)
        if key not in remote or version(filename) > version(remote[key][0]):
            remote[key] = (filename, address)

    if not remote and not known:
        raise FileNotFoundError(f"No full-cadence MAG RTN files listed for {day}")

    # Only download selected highest-version segments that we don't have.
    for key, (filename, address) in sorted(remote.items()):
        if key in known and version(known[key].name) >= version(filename):
            continue
        target = FULL_CACHE / filename
        partial = target.with_suffix(".cdf.part")
        try:
            with requests.get(address, stream=True, timeout=(30, 180)) as stream:
                stream.raise_for_status()
                with partial.open("wb") as out:
                    for block in stream.iter_content(1024 * 1024):
                        if block:
                            out.write(block)
            if partial.stat().st_size < 1024:
                raise ValueError(f"Download unexpectedly small: {filename}")
            partial.replace(target)
        finally:
            partial.unlink(missing_ok=True)
        known[key] = target

    # If the NASA inventory lists fewer than four segments, scanning the
    # available coverage is allowed, but missing coverage must be disclosed.
    paths = [known[key] for key in sorted(known)]
    missing = sorted(expected - set(known))
    if missing:
        print(f"[FIELDS full cadence] {day}: missing segments {missing}; "
              "screening only available coverage")
    return paths


def iter_full_cadence(day, records_per_chunk=250_000):
    """Yield chronological, quality-screened full-cadence magnetic DataFrames.

    Uses CDF record ranges to bound RAM. The CDF files themselves are cached,
    never deleted based on a candidate-screening decision.
    """
    if records_per_chunk < 1000:
        raise ValueError("records_per_chunk must be at least 1000")
    for path in _full_day_files(day):
        cdf = cdflib.CDF(str(path))
        names = {name.lower(): name for name in (
            list(cdf.cdf_info().zVariables) + list(cdf.cdf_info().rVariables)
        )}
        if FULL_VECTOR.lower() not in names or FULL_EPOCH.lower() not in names:
            raise ValueError(f"Missing full-cadence vector/epoch in {path.name}")
        vector = names[FULL_VECTOR.lower()]
        epoch = names[FULL_EPOCH.lower()]
        last = cdf.varinq(vector).Last_Rec
        fill = cdf.varattsget(vector).get("FILLVAL")

        flags = None
        if {"psp_fld_l2_quality_flags", "epoch_quality_flags"}.issubset(names):
            flag_dates = pd.to_datetime(cdflib.cdfepoch.to_datetime(
                cdf.varget(names["epoch_quality_flags"])
            ))
            flag_values = np.asarray(
                cdf.varget(names["psp_fld_l2_quality_flags"])
            ).ravel()
            if len(flag_dates) != len(flag_values):
                raise ValueError(f"Invalid quality flags in {path.name}")
            flags = pd.Series(flag_values, index=flag_dates).sort_index()
            flags = flags[~flags.index.duplicated(keep="last")]

        for first in range(0, last + 1, records_per_chunk):
            final = min(first + records_per_chunk - 1, last)
            raw = np.asarray(
                cdf.varget(vector, startrec=first, endrec=final), dtype=float
            )
            index = pd.to_datetime(cdflib.cdfepoch.to_datetime(
                cdf.varget(epoch, startrec=first, endrec=final)
            ))
            if raw.ndim != 2 or raw.shape[1] != 3 or len(raw) != len(index):
                raise ValueError(f"Unexpected full-cadence vector shape in {path.name}")
            if fill is not None and np.isscalar(fill):
                raw[raw == float(fill)] = np.nan
            raw[(~np.isfinite(raw)) | (np.abs(raw) > 1e6)] = np.nan
            frame = pd.DataFrame(raw, index=index, columns=["BR", "BT", "BN"])
            frame = frame[~frame.index.duplicated(keep="last")].sort_index()
            if flags is not None and not frame.empty:
                aligned = flags.reindex(
                    frame.index, method="nearest",
                    tolerance=pd.Timedelta(seconds=35)
                )
                bad = aligned.isna().to_numpy() | (
                    (aligned.fillna(0).astype("int64").to_numpy() & BAD_BITS) != 0
                )
                frame.loc[bad, ["BR", "BT", "BN"]] = np.nan
            frame["Bmag"] = np.sqrt(
                (frame[["BR", "BT", "BN"]] ** 2).sum(axis=1, min_count=3)
            )
            frame.index.name = "UTC"
            if not frame.empty:
                yield frame


def full_cadence_seconds(day, records_per_chunk=250_000):
    """Stream raw samples into a 1-second *screening* representation.

    Carry the unfinished final second into the next chunk to avoid
    sample-bin artifacts at arbitrary record boundaries.
    """
    seconds = []
    pending = None
    for chunk in iter_full_cadence(day, records_per_chunk):
        frame = (pd.concat([pending, chunk]) if pending is not None else chunk)
        end_second = frame.index[-1].floor("s")
        done = frame.loc[frame.index < end_second, ["BR", "BT", "BN"]]
        pending = frame.loc[frame.index >= end_second, ["BR", "BT", "BN"]]
        if not done.empty:
            seconds.append(done.resample("1s").median())
    if pending is not None and not pending.empty:
        seconds.append(pending.resample("1s").median())
    if not seconds:
        raise RuntimeError(f"No full-cadence FIELDS measurements for {day}")
    result = pd.concat(seconds).groupby(level=0).median().sort_index()
    result["Bmag"] = np.sqrt(
        (result[["BR", "BT", "BN"]] ** 2).sum(axis=1, min_count=3)
    )
    result.index.name = "UTC"
    return result


def pre_screen_magnetic_rotations(
    seconds, angle_threshold=40.0, peak_threshold=50.0,
    minimum_seconds=3, context_minutes=5
):
    """Permissive magnetic-rotation pre-screen; NOT event classification.

    Local direction uses a centered 181-minute median of *minute summaries
    derived from full-cadence samples*. Candidate windows retain surrounding
    context for later detailed physical interpretation.
    """
    empty = ["start", "end", "category", "peak_angle_deg",
             "core_start", "core_end", "screening_method"]
    if seconds.empty:
        return pd.DataFrame(columns=empty)
    b = seconds[["BR", "BT", "BN"]].sort_index()
    b = b[~b.index.duplicated(keep="last")]
    if len(b) < 120:
        return pd.DataFrame(columns=empty)
    # Avoid reference collapse due to transient vector reversals. No
    # normalization is attempted until a nonzero reference is available.
    minute = b.resample("1min").median()
    bg = minute.rolling("181min", center=True, min_periods=30).median()
    ref = bg.reindex(b.index, method="nearest",
                     tolerance=pd.Timedelta(seconds=35))
    vectors = b.to_numpy(dtype=float)
    background = ref.to_numpy(dtype=float)
    denominator = (np.linalg.norm(vectors, axis=1) *
                   np.linalg.norm(background, axis=1))
    ok = (np.isfinite(vectors).all(axis=1) &
          np.isfinite(background).all(axis=1) &
          np.isfinite(denominator) & (denominator > 1e-8))
    cosine = np.full(len(b), np.nan)
    cosine[ok] = (
        np.sum(vectors[ok] * background[ok], axis=1) / denominator[ok]
    )
    angle = pd.Series(
        np.degrees(np.arccos(np.clip(cosine, -1, 1))), index=b.index
    )
    # Three seconds above the permissive threshold within a five-second
    # window. This suppresses isolated bad readings, without the strict
    # duration rule used for the later event analysis.
    above = angle.ge(angle_threshold)
    supported = (above.astype(int).rolling(
        5, center=True, min_periods=1
    ).sum() >= minimum_seconds) & angle.notna()
    positions = np.flatnonzero(supported.to_numpy())
    if not len(positions):
        return pd.DataFrame(columns=empty)
    groups = np.split(positions, np.where(
        np.diff(positions) > 1
    )[0] + 1)
    spans = [(int(group[0]), int(group[-1])) for group in groups]
    merged = [spans[0]]
    for start, end in spans[1:]:
        previous_start, previous_end = merged[-1]
        dt = (b.index[start] - b.index[previous_end]).total_seconds()
        # Don't bridge real missing-data intervals.
        if dt <= 20 and angle.iloc[previous_end:start+1].notna().all():
            merged[-1] = (previous_start, end)
        else:
            merged.append((start, end))
    padding = pd.Timedelta(minutes=context_minutes)
    results = []
    for first, last in merged:
        core = angle.iloc[first:last + 1]
        if int(above.iloc[first:last+1].sum()) < minimum_seconds:
            continue
        peak = core.max()
        if not np.isfinite(peak) or peak < peak_threshold:
            continue
        left = max(b.index[0], b.index[first] - padding)
        right = min(b.index[-1], b.index[last] + padding)
        row = dict(
            start=left, end=right, category="screened_full",
            peak_angle_deg=round(float(peak), 1),
            core_start=b.index[first], core_end=b.index[last],
            screening_method="full_cadence_1s_181min_local"
        )
        if results and left <= results[-1]["end"]:
            results[-1]["end"] = max(results[-1]["end"], right)
            results[-1]["core_end"] = row["core_end"]
            results[-1]["peak_angle_deg"] = max(
                results[-1]["peak_angle_deg"], row["peak_angle_deg"]
            )
        else:
            results.append(row)
    return pd.DataFrame(results, columns=empty)


def screen_full_cadence(day, records_per_chunk=250_000):
    """Scan all available full-cadence segments and return candidate windows."""
    seconds = full_cadence_seconds(day, records_per_chunk)
    return pre_screen_magnetic_rotations(seconds)


def read_full_interval(start, end, records_per_chunk=250_000):
    """Yield only original full-cadence samples inside an observation window.

    Call this from the detailed analysis stage; no interpolation is performed.
    """
    left, right = pd.Timestamp(start), pd.Timestamp(end)
    if left.tzinfo is not None:
        left = left.tz_convert("UTC").tz_localize(None)
    if right.tzinfo is not None:
        right = right.tz_convert("UTC").tz_localize(None)
    if right < left:
        raise ValueError("end must not precede start")
    for day in pd.date_range(left.normalize(), right.normalize()):
        for chunk in iter_full_cadence(day.date(), records_per_chunk):
            selected = chunk.loc[(chunk.index >= left) & (chunk.index <= right)]
            if not selected.empty:
                yield selected


def _screen_self_test():
    """Synthetic exercise of the reader's permissive pre-screening logic."""
    idx = pd.date_range("2018-11-06", periods=10800, freq="s")
    vectors = pd.DataFrame(
        {"BR": np.full(len(idx), 80.0),
         "BT": np.zeros(len(idx)), "BN": np.zeros(len(idx))},
        index=idx
    )
    assert pre_screen_magnetic_rotations(vectors).empty
    vectors.loc["2018-11-06 01:30:00":"2018-11-06 01:30:59", "BR"] = -80.0
    found = pre_screen_magnetic_rotations(vectors)
    assert len(found) == 1, f"Expected one synthetic rotation: {found}"
    assert found.iloc[0]["core_start"] <= pd.Timestamp("2018-11-06 01:30:05")
    assert found.iloc[0]["core_end"] >= pd.Timestamp("2018-11-06 01:30:55")
    print("PASS: full-cadence pre-screen synthetic quiet/rotation cases")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Read and pre-screen PSP/FIELDS RTN data")
    parser.add_argument("--start", default="2018-11-06")
    parser.add_argument("--days", type=int, default=1)
    parser.add_argument("--resolution", choices=tuple(PRODUCTS), default="1min",
                        help="Existing direct readers; full cadence is streamed via --screen-full")
    parser.add_argument("--screen-full", action="store_true",
                        help="Stream full cadence, perform permissive pre-screen, save candidate windows")
    parser.add_argument("--self-test", action="store_true",
                        help="Run synthetic pre-screen test without downloading NASA data")
    args = parser.parse_args()
    if args.self_test:
        _screen_self_test()
    else:
        for day in days(args.start, args.days):
            if args.screen_full:
                candidates = screen_full_cadence(day)
                output = Path("data") / f"full_screen_{day:%Y%m%d}.csv"
                output.parent.mkdir(parents=True, exist_ok=True)
                candidates.to_csv(output, index=False)
                print(f"{day}: {len(candidates)} pre-screened magnetic windows; "
                      f"saved {output}")
            else:
                df = get_day(day, resolution=args.resolution)
                print(f"{day} [{args.resolution}]: {df.Bmag.notna().sum()} "
                      f"valid samples of {len(df)}")
