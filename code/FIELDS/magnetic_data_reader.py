"""Parker Solar Probe FIELDS magnetic-field reader.

One interface for PSP/FIELDS Level-2 RTN products:
  get_day(date, resolution="1min")   -- one-minute averages
  get_day(date, resolution="4cyc")   -- four samples per instrument cycle
  scan_time_range(start, end)        -- arbitrary interval, full-cadence pre-screen
  read_full_interval(start, end)     -- original full-cadence samples

The direct readers return UTC-indexed BR, BT, BN and Bmag (nT).
The products are not interchangeable: the low-cadence product is not
interpolated into high-cadence measurements. The reader only performs a
permissive pre-screen. Scientific characterization belongs to
magnetic_data_analysis.py.
"""

import argparse
import json
from tempfile import TemporaryDirectory
from datetime import date, timedelta
from pathlib import Path
from urllib.parse import urlparse
import re
import sys
import time

import cdflib
import numpy as np
import pandas as pd
import requests

BASE_API = "https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/sp_phys/datasets"

def _progress_text(label, current, total=None, unit="records"):
    """Compact ASCII bar for known totals; byte counter if total is unknown."""
    if total is not None and total > 0:
        fraction = min(1.0, max(0.0, current / total))
        filled = int(28 * fraction)
        track = "#" * filled + "-" * (28 - filled)
        if unit == "bytes":
            amounts = f"{current / 1048576:.1f}/{total / 1048576:.1f} MiB"
        else:
            amounts = f"{current:,}/{total:,} {unit}"
        return f"[{label}] [{track}] {fraction * 100:5.1f}%  {amounts}"
    amount = (f"{current / 1048576:.1f} MiB received" if unit == "bytes"
              else f"{current:,} {unit}")
    return f"[{label}] {amount} (total size unknown)"


class _ProgressBar:
    """Live stderr progress without any extra package requirement."""

    def __init__(self, label, total=None, unit="records"):
        self.label = label
        self.total = total
        self.unit = unit
        self.current = 0
        self.last_draw = 0.0
        self.last_width = 0
        self.draw(force=True)

    def draw(self, force=False):
        now = time.monotonic()
        interactive = sys.stderr.isatty()
        # A live terminal gets an in-place bar; redirected/background logs
        # get readable progress lines about once a second.
        interval = 0.2 if interactive else 1.0
        if not force and now - self.last_draw < interval:
            return
        message = _progress_text(
            self.label, self.current, self.total, self.unit
        )
        if interactive:
            sys.stderr.write("\r" + message.ljust(self.last_width))
            self.last_width = max(self.last_width, len(message))
        else:
            sys.stderr.write(message + "\n")
        sys.stderr.flush()
        self.last_draw = now

    def advance(self, count):
        self.current += count
        self.draw()

    def close(self):
        self.draw(force=True)
        if sys.stderr.isatty():
            sys.stderr.write("\n")
            sys.stderr.flush()


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


def _as_utc_naive(value):
    """Normalize an ISO date/time to UTC for matching NASA UTC CDF epochs."""
    timestamp = pd.Timestamp(value)
    if pd.isna(timestamp):
        raise ValueError("A valid timestamp is required")
    return (timestamp.tz_convert("UTC").tz_localize(None)
            if timestamp.tzinfo is not None else timestamp)


def _segment_hours(day, start=None, end=None):
    """Return six-hour files overlapping the requested [start, end) interval."""
    midnight = pd.Timestamp(day).normalize()
    left = _as_utc_naive(start) if start is not None else midnight
    right = _as_utc_naive(end) if end is not None else midnight + pd.Timedelta(days=1)
    return tuple(hour for hour in (0, 6, 12, 18)
                 if midnight + pd.Timedelta(hours=hour) < right
                 and midnight + pd.Timedelta(hours=hour+6) > left)




def _nasa_catalog_json(url, attempts=6, pause=time.sleep):
    """Retry transient CDAWeb metadata API failures as well as CDF transfers."""
    for attempt in range(1, attempts + 1):
        try:
            response = requests.get(
                url, headers={"Accept": "application/json"}, timeout=90)
            response.raise_for_status()
            result = response.json()
            if not isinstance(result, dict):
                raise ValueError("NASA response is not a JSON object")
            return result
        except (requests.exceptions.RequestException, ValueError) as exc:
            status = (exc.response.status_code
                      if isinstance(exc, requests.exceptions.HTTPError)
                      and exc.response is not None else None)
            if (status is not None and status not in (408, 429)
                    and status < 500):
                raise
            if attempt == attempts:
                raise RuntimeError(
                    f"NASA catalog unavailable after {attempts} attempts: "
                    f"{exc}") from exc
            wait = min(45, 2 ** (attempt - 1) * 3)
            print(f"[NASA RETRY] Metadata request {attempt}/{attempts}: "
                  f"{exc}; retrying in {wait}s", flush=True)
            pause(wait)
    raise RuntimeError("Unreachable NASA catalog retry state")


def _download_full_cdf(address, target, attempts=6, pause=time.sleep):
    """Retry interrupted NASA CDF downloads and resume verified byte ranges.

    A .cdf.part file is intentionally kept on transient failures. We publish
    the finished .cdf only after verifying the transfer size and CDF structure.
    If the server ignores Range, restart from zero rather than append a
    duplicate file. All failures remain visible in the scan log.
    """
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(".cdf.part")
    for attempt in range(1, attempts + 1):
        existing = partial.stat().st_size if partial.exists() else 0
        headers = {"Accept-Encoding": "identity"}
        if existing:
            headers["Range"] = f"bytes={existing}-"
        try:
            print(f"[NASA] Downloading {target.name} "
                  f"(attempt {attempt}/{attempts}, "
                  f"{existing/1048576:.1f} MiB already cached)", flush=True)
            with requests.get(address, headers=headers, stream=True,
                              timeout=(30, 180)) as stream:
                if stream.status_code == 416 and existing:
                    # Cached partial is larger than the server's object.
                    partial.unlink(missing_ok=True)
                    raise requests.exceptions.ConnectionError(
                        "Cached partial exceeds remote file; restarting")
                stream.raise_for_status()

                use_append = existing > 0 and stream.status_code == 206
                expected_total = None
                if use_append:
                    content_range = stream.headers.get("Content-Range", "")
                    match = re.fullmatch(
                        r"bytes (\d+)-(\d+)/(\d+)", content_range)
                    if (match is None or int(match.group(1)) != existing
                            or int(match.group(2)) < existing):
                        partial.unlink(missing_ok=True)
                        raise requests.exceptions.ConnectionError(
                            "Invalid resumed Content-Range; restarting")
                    expected_total = int(match.group(3))
                    if expected_total <= existing:
                        partial.unlink(missing_ok=True)
                        raise requests.exceptions.ConnectionError(
                            "Invalid resumed object length; restarting")
                    mode = "ab"
                elif stream.status_code == 200:
                    # The server did not honor Range: overwrite, never append.
                    mode = "wb"
                else:
                    raise requests.exceptions.ConnectionError(
                        f"Unexpected download status: {stream.status_code}")

                length = stream.headers.get("Content-Length")
                if length is not None and length.isdigit():
                    received_size = int(length)
                    total_from_length = (
                        existing + received_size if use_append
                        else received_size)
                    if expected_total is not None and (
                            total_from_length != expected_total):
                        partial.unlink(missing_ok=True)
                        raise requests.exceptions.ConnectionError(
                            "Inconsistent response length for resumed CDF")
                    expected_total = total_from_length

                progress = _ProgressBar(
                    "Download", total=expected_total, unit="bytes")
                if mode == "ab":
                    progress.advance(existing)
                try:
                    with partial.open(mode) as out:
                        for block in stream.iter_content(1024 * 1024):
                            if block:
                                out.write(block)
                                progress.advance(len(block))
                finally:
                    progress.close()

            if not partial.is_file() or partial.stat().st_size < 1024:
                raise IOError(f"Downloaded CDF too small: {target.name}")
            if expected_total is not None and (
                    partial.stat().st_size != expected_total):
                raise requests.exceptions.ConnectionError(
                    f"Incomplete CDF: got {partial.stat().st_size}, "
                    f"expected {expected_total} bytes")

            # Check that NASA's CDF structure is actually readable before
            # promoting a partial download into the permanent cache.
            try:
                cdf = cdflib.CDF(str(partial))
                info = cdf.cdf_info()
                variables = {v.lower() for v in (
                    list(info.zVariables) + list(info.rVariables))}
                if (FULL_VECTOR.lower() not in variables or
                        FULL_EPOCH.lower() not in variables):
                    raise ValueError("Expected FIELDS variables missing")
            except Exception as exc:
                partial.unlink(missing_ok=True)
                raise ValueError(
                    f"Downloaded object is not a valid FIELDS CDF: {exc}"
                ) from exc

            partial.replace(target)
            print(f"[NASA] Saved {target.name}", flush=True)
            return target
        except (requests.exceptions.RequestException, OSError,
                ValueError) as exc:
            if attempt == attempts:
                raise RuntimeError(
                    f"Failed to download {target.name} after {attempts} "
                    f"attempts: {exc}") from exc
            wait = min(45, 2 ** (attempt - 1) * 3)
            print(f"[NASA RETRY] {target.name}: {type(exc).__name__}: "
                  f"{exc}; retrying in {wait}s", flush=True)
            pause(wait)
    raise RuntimeError(f"Unreachable download state: {target.name}")


def _full_day_files(day, start=None, end=None):
    """Cache only full-cadence six-hour segments intersecting [start, end)."""
    day = pd.Timestamp(day).date()
    marker = day.strftime("%Y%m%d")
    required = {f"{marker}{hour:02d}" for hour in
                _segment_hours(day, start=start, end=end)}
    if not required:
        return []
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

    if required.issubset(known):
        return [known[key] for key in sorted(required)]

    stop = day + timedelta(days=1)
    api = f"{BASE_API}/{FULL_DATASET}/orig_data"
    url = f"{api}/{day:%Y%m%d}T000000Z,{stop:%Y%m%d}T000000Z/"
    print(f"[NASA] Looking up full-cadence files for {day}...", flush=True)
    descriptions = _nasa_catalog_json(url).get("FileDescription", [])
    if isinstance(descriptions, dict):
        descriptions = [descriptions]
    remote = {}
    for entry in descriptions:
        address = entry.get("Name", "")
        filename = Path(urlparse(address).path).name
        match = FULL_NAME.fullmatch(filename)
        if not match or match.group(1) not in required:
            continue
        key = match.group(1)
        if key not in remote or version(filename) > version(remote[key][0]):
            remote[key] = (filename, address)

    if not remote and not (required & set(known)):
        raise FileNotFoundError(f"No requested full-cadence MAG RTN files listed for {day}")

    # Only download selected highest-version segments that we don't have.
    for key, (filename, address) in sorted(remote.items()):
        if key in known and version(known[key].name) >= version(filename):
            continue
        target = FULL_CACHE / filename
        known[key] = _download_full_cdf(address, target)

    # Requested 6-hour files may be absent; any incomplete scan coverage
    # is explicitly reported to avoid treating data gaps as quiet solar wind.
    paths = [known[key] for key in sorted(required & set(known))]
    missing = sorted(required - set(known))
    if missing:
        print(f"[FIELDS full cadence] {day}: missing segments {missing}; "
              "screening only available coverage")
    return paths


def _lower_record_bound(last_record, timestamp_at, target):
    """First record at/after target; efficient for short selections in huge CDFs.

    NASA documents MAG epoch timestamps as monotonically increasing.
    Only O(log N) timestamp records need reading, instead of decoding all
    magnetic vectors in a six-hour file to inspect a few requested minutes.
    """
    low, high = 0, last_record + 1
    while low < high:
        mid = (low + high) // 2
        if timestamp_at(mid) < target:
            low = mid + 1
        else:
            high = mid
    return low


def iter_full_cadence(day, records_per_chunk=250_000, start=None, end=None):
    """Yield chronological, quality-screened full-cadence magnetic DataFrames.

    Uses CDF record ranges to bound RAM. The CDF files themselves are cached,
    never deleted based on a candidate-screening decision.
    """
    if records_per_chunk < 1000:
        raise ValueError("records_per_chunk must be at least 1000")
    left = _as_utc_naive(start) if start is not None else None
    right = _as_utc_naive(end) if end is not None else None
    for path in _full_day_files(day, start=left, end=right):
        cdf = cdflib.CDF(str(path))
        names = {name.lower(): name for name in (
            list(cdf.cdf_info().zVariables) + list(cdf.cdf_info().rVariables)
        )}
        if FULL_VECTOR.lower() not in names or FULL_EPOCH.lower() not in names:
            raise ValueError(f"Missing full-cadence vector/epoch in {path.name}")
        vector = names[FULL_VECTOR.lower()]
        epoch = names[FULL_EPOCH.lower()]
        last = int(cdf.varinq(vector).Last_Rec)
        if last < 0:
            continue
        fill = cdf.varattsget(vector).get("FILLVAL")

        # Source filenames identify 6-hour UTC intervals. Avoid reading
        # vectors before/after a short requested interval. For scans that
        # cover the whole file we simply stream all records.
        filename_match = FULL_NAME.fullmatch(path.name)
        if filename_match is None:
            raise ValueError(f"Unrecognized full-cadence filename: {path.name}")
        segment_start = pd.to_datetime(filename_match.group(1), format="%Y%m%d%H")
        segment_end = segment_start + pd.Timedelta(hours=6)

        def timestamp_at(record):
            one = cdf.varget(epoch, startrec=record, endrec=record)
            converted = cdflib.cdfepoch.to_datetime(one)
            return pd.Timestamp(np.atleast_1d(converted)[0])

        first_record = (0 if left is None or left <= segment_start else
                        _lower_record_bound(last, timestamp_at, left))
        stop_record = (last + 1 if right is None or right >= segment_end else
                       _lower_record_bound(last, timestamp_at, right))
        if first_record >= stop_record:
            continue

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

        print(f"[FIELDS] Reading {path.name}", flush=True)
        progress = _ProgressBar("Processing", total=stop_record-first_record)
        try:
            for first in range(first_record, stop_record, records_per_chunk):
                final = min(first + records_per_chunk - 1, stop_record - 1)
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
                if left is not None:
                    frame = frame.loc[frame.index >= left]
                if right is not None:
                    frame = frame.loc[frame.index < right]
                if not frame.empty:
                    if flags is not None:
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
                progress.advance(final-first+1)
                if not frame.empty:
                    yield frame
        finally:
            progress.close()


def full_cadence_seconds(day, records_per_chunk=250_000, start=None, end=None):
    """Stream raw samples into a 1-second *screening* representation.

    Carry the unfinished final second into the next chunk to avoid
    sample-bin artifacts at arbitrary record boundaries.
    """
    seconds = []
    pending = None
    for chunk in iter_full_cadence(day, records_per_chunk, start=start, end=end):
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
    context for later detailed physical interpretation. The 1-second medians
    can miss sub-second changes; this is deliberately an exploratory screen,
    not an exhaustive scientific switchback detector.
    """
    empty = ["start", "end", "category", "peak_angle_deg",
             "core_start", "core_end", "screening_method"]
    if seconds.empty:
        return pd.DataFrame(columns=empty)
    b = seconds[["BR", "BT", "BN"]].sort_index()
    b = b[~b.index.duplicated(keep="last")]
    # A missing telemetry interval must be represented by NaNs. Without
    # explicit one-second bins, two widely separated valid samples could
    # masquerade as a continuous event and cross a real data gap.
    b = b.resample("1s").median()
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
        for chunk in iter_full_cadence(day.date(), records_per_chunk, start=left, end=right):
            selected = chunk.loc[(chunk.index >= left) & (chunk.index <= right)]
            if not selected.empty:
                yield selected



def _parse_duration(start, duration):
    """A single request may specify minutes, hours, days, weeks or months."""
    match = re.fullmatch(r"(\d+)(mo|min|m|h|d|w)", duration.strip().lower())
    if not match:
        raise ValueError("Duration must look like 10m, 2h, 5d, 3w or 2mo")
    count, unit = int(match.group(1)), match.group(2)
    if count == 0:
        raise ValueError("Duration must be positive")
    if unit == "mo":
        return start + pd.DateOffset(months=count)
    return start + pd.Timedelta(**{
        "m": {"minutes": count}, "min": {"minutes": count},
        "h": {"hours": count}, "d": {"days": count}, "w": {"weeks": count}
    }[unit])


def _merge_screened_windows(records):
    """Combine overlapping context windows; event times are provisional."""
    if not records:
        return pd.DataFrame(columns=[
            "start", "end", "category", "peak_angle_deg",
            "core_start", "core_end", "screening_method"
        ])
    rows = sorted(records, key=lambda row: row["start"])
    result = [rows[0].copy()]
    for row in rows[1:]:
        previous = result[-1]
        if row["start"] <= previous["end"]:
            previous["end"] = max(previous["end"], row["end"])
            previous["core_start"] = min(previous["core_start"], row["core_start"])
            previous["core_end"] = max(previous["core_end"], row["core_end"])
            previous["peak_angle_deg"] = max(
                previous["peak_angle_deg"], row["peak_angle_deg"])
        else:
            result.append(row.copy())
    return pd.DataFrame(result)



# Incremental catalogs: one CSV per user-requested range, updated after
# every completed day. An internal JSON checkpoint supports interruption.
SCAN_STATE_VERSION = 1
CANDIDATE_COLUMNS = [
    "start", "end", "category", "peak_angle_deg",
    "core_start", "core_end", "screening_method"
]
TIME_COLUMNS = ("start", "end", "core_start", "core_end")


def _catalog_from_records(records):
    merged = _merge_screened_windows(records)
    if merged.empty:
        return pd.DataFrame(columns=CANDIDATE_COLUMNS)
    return merged.reindex(columns=CANDIDATE_COLUMNS).sort_values(
        "start"
    ).reset_index(drop=True)


def _atomic_text(path, content):
    """Write atomically so VS Code never opens a half-written catalog."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            handle.write(content)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _checkpoint_file(output):
    return Path(str(output) + ".checkpoint.json")


def _save_scan_checkpoint(output, start, end, completed_day, rows, valid_seconds):
    """Save restart state first, then publish the candidate CSV.

    If interrupted after checkpointing but before the CSV update,
    --resume repairs the CSV automatically from checkpoint state.
    """
    output = Path(output)
    state = {
        "version": SCAN_STATE_VERSION,
        "start": _as_utc_naive(start).isoformat(),
        "end": _as_utc_naive(end).isoformat(),
        "completed_day": pd.Timestamp(completed_day).date().isoformat(),
        "valid_seconds": int(valid_seconds),
        "raw_candidates": [
            {key: (pd.Timestamp(value).isoformat() if key in TIME_COLUMNS
                   else value) for key, value in row.items()}
            for row in rows
        ],
    }
    _atomic_text(_checkpoint_file(output), json.dumps(state, indent=2) + "\n")
    catalog = _catalog_from_records(rows)
    _atomic_text(output, catalog.to_csv(index=False))
    print(f"[SAVED] {state['completed_day']}: {len(catalog)} windows "
          f"available now in {output}", flush=True)
    return catalog


def _load_scan_checkpoint(output, start, end):
    """Load saved progress for precisely the same requested UTC interval."""
    state_path = _checkpoint_file(output)
    if not state_path.is_file():
        raise FileNotFoundError(
            f"No checkpoint at {state_path}. Run without --resume to start. "
            "Scans started before checkpointing was added cannot be resumed."
        )
    state = json.loads(state_path.read_text(encoding="utf-8"))
    if (state.get("version") != SCAN_STATE_VERSION
        or state.get("start") != _as_utc_naive(start).isoformat()
        or state.get("end") != _as_utc_naive(end).isoformat()):
        raise ValueError(
            "Checkpoint does not match this scan or code version. "
            "Use the original start/end/output or start without --resume."
        )
    completed_day = pd.Timestamp(state["completed_day"]).normalize()
    last_day = (_as_utc_naive(end) - pd.Timedelta(nanoseconds=1)).normalize()
    if not _as_utc_naive(start).normalize() <= completed_day <= last_day:
        raise ValueError("Checkpoint day falls outside the requested interval")
    rows = []
    for saved in state["raw_candidates"]:
        row = dict(saved)
        for key in TIME_COLUMNS:
            row[key] = _as_utc_naive(row[key])
        rows.append(row)
    valid_seconds = int(state["valid_seconds"])
    if valid_seconds < 0:
        raise ValueError("Invalid checkpoint coverage")
    # Repair stale CSV after a crash between the two atomic writes.
    _atomic_text(output, _catalog_from_records(rows).to_csv(index=False))
    return completed_day, rows, valid_seconds


def scan_time_range(start, end, records_per_chunk=250_000, on_day=None,
                    on_checkpoint=None, resume_after=None,
                    initial_records=None, initial_valid_seconds=0):
    """Screen an arbitrary UTC [start, end) interval at full MAG cadence.

    The user specifies one time range; the implementation automatically
    processes each day with 100-minute context across date boundaries. Only
    one day's 1-second summaries plus neighboring margins are in RAM at a time.
    Every day is read from original full-cadence CDFs, never from the NASA
    1-minute or 4-samples/cycle products.

    Output windows are candidate *investigation* intervals, not physically
    validated switchback boundaries. Missing segments and sparse coverage are
    reported; scans continue across unavailable days and fail if no usable
    magnetic samples exist anywhere in the requested interval.
    """
    left, right = _as_utc_naive(start), _as_utc_naive(end)
    if right <= left:
        raise ValueError("End time must be later than start time")
    halo = pd.Timedelta(minutes=100)
    rows = list(initial_records) if initial_records is not None else []
    total_valid_seconds = int(initial_valid_seconds)
    completed_before = (_as_utc_naive(resume_after).normalize()
                        if resume_after is not None else None)
    last_day = (right - pd.Timedelta(nanoseconds=1)).normalize()
    all_days = pd.date_range(left.normalize(), last_day, freq="D")
    print(f"[SCAN] {left} to {right} UTC ({len(all_days)} calendar day(s))",
          flush=True)
    for day_number, day in enumerate(all_days, start=1):
        if completed_before is not None and day <= completed_before:
            continue
        print(f"[SCAN] Day {day_number}/{len(all_days)}: {day.date()}",
              flush=True)
        core_left = max(left, day)
        core_right = min(right, day + pd.Timedelta(days=1))
        expected_seconds = max(1, int(np.ceil((core_right - core_left).total_seconds())))
        around_left, around_right = core_left - halo, core_right + halo
        parts = []
        halo_last = (around_right - pd.Timedelta(nanoseconds=1)).normalize()
        for neighboring_day in pd.date_range(
            around_left.normalize(), halo_last, freq="D"
        ):
            midnight = neighboring_day
            piece_start = max(around_left, midnight)
            piece_end = min(around_right, midnight + pd.Timedelta(days=1))
            try:
                part = full_cadence_seconds(
                    midnight.date(), records_per_chunk,
                    start=piece_start, end=piece_end
                )
            except FileNotFoundError:
                # Coverage can be intermittent. Record the gap, continue to
                # the next available segment, and NEVER label a gap quiet.
                print(f"[MISSING COVERAGE] {piece_start} to {piece_end}: "
                      "NASA full-cadence file not available", flush=True)
                continue
            except RuntimeError as exc:
                if "No full-cadence FIELDS measurements" not in str(exc):
                    raise
                print(f"[MISSING COVERAGE] {piece_start} to {piece_end}: "
                      "no measurements in available files", flush=True)
                continue
            parts.append(part)
        if not parts:
            print(_progress_text("Overall scan", day_number, len(all_days), "days"),
                  flush=True)
            print(f"[COVERAGE WARNING] {core_left} to {core_right}: no valid "
                  "observations (not a quiet interval)", flush=True)
            if on_day is not None:
                on_day(day.date(), 0, 0, expected_seconds)
            if on_checkpoint is not None:
                on_checkpoint(day.date(), rows, total_valid_seconds)
            continue
        series = pd.concat(parts).sort_index()
        series = series[~series.index.duplicated(keep="last")]
        core = series.loc[(series.index >= core_left) & (series.index < core_right)]
        if core.empty or not core[["BR","BT","BN"]].notna().all(axis=1).any():
            print(_progress_text("Overall scan", day_number, len(all_days), "days"),
                  flush=True)
            print(f"[COVERAGE WARNING] {core_left} to {core_right}: no valid "
                  "magnetic vectors (not a quiet interval)", flush=True)
            if on_day is not None:
                on_day(day.date(), 0, 0, expected_seconds)
            if on_checkpoint is not None:
                on_checkpoint(day.date(), rows, total_valid_seconds)
            continue
        print("[SCREEN] Checking magnetic rotations...", flush=True)
        candidate = pre_screen_magnetic_rotations(series)
        if not candidate.empty:
            # Prevent duplication from the 100-minute context on adjacent days.
            candidate = candidate.loc[
                (candidate["core_end"] >= core_left) &
                (candidate["core_start"] < core_right)
            ]
            for row in candidate.to_dict("records"):
                rows.append(row)
        observed = int(core[["BR","BT","BN"]].notna().all(axis=1).sum())
        total_valid_seconds += observed
        coverage = min(1.0, observed / expected_seconds)
        if coverage < 0.90:
            print(f"[COVERAGE WARNING] {day.date()}: {coverage:.1%} of requested "
                  "one-second bins have valid magnetic vectors; missing time "
                  "must not be interpreted as an absence of switchbacks")
        if on_day is not None:
            on_day(day.date(), len(candidate), observed, expected_seconds)
        print(_progress_text("Overall scan", day_number, len(all_days), "days"),
              flush=True)
        if on_checkpoint is not None:
            on_checkpoint(day.date(), rows, total_valid_seconds)
    if total_valid_seconds == 0:
        raise RuntimeError("No valid full-cadence observations in the requested "
                           "interval. No switchback screening catalog produced.")
    return _catalog_from_records(rows)


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
    assert _segment_hours("2018-11-06",
                          "2018-11-06T10:12:00",
                          "2018-11-06T10:22:00") == (6,)
    assert _segment_hours("2018-11-06",
                          "2018-11-06T05:55:00",
                          "2018-11-06T06:05:00") == (0, 6)
    t0 = pd.Timestamp("2018-11-06T00:00:00")
    fake_epoch = [t0 + pd.Timedelta(milliseconds=100*i) for i in range(100)]
    lookup = lambda i: fake_epoch[i]
    assert _lower_record_bound(99, lookup, t0) == 0
    assert _lower_record_bound(99, lookup, t0 + pd.Timedelta(milliseconds=250)) == 3
    assert _lower_record_bound(99, lookup, t0 + pd.Timedelta(seconds=10)) == 100
    # Sparse observations should not bridge missing-data intervals.
    gap_index = pd.date_range("2018-11-06", periods=10800, freq="s")
    sparse = pd.DataFrame(
        {"BR": 80., "BT": 0., "BN": 0.}, index=gap_index
    )
    sparse.loc["2018-11-06 01:30:00":"2018-11-06 01:30:09", "BR"] = -80.
    sparse.loc["2018-11-06 01:32:00":"2018-11-06 01:32:09", "BR"] = -80.
    sparse = sparse.drop(index=pd.date_range(
        "2018-11-06 01:30:10", "2018-11-06 01:31:59", freq="s"))
    observed = pre_screen_magnetic_rotations(sparse, context_minutes=0)
    assert len(observed) == 2, f"Telemetric data gaps were bridged: {observed}"
    assert _parse_duration(pd.Timestamp("2018-11-06T12:00:00"), "10m") == pd.Timestamp("2018-11-06T12:10:00")
    assert _parse_duration(pd.Timestamp("2018-11-06T12:00:00"), "2mo") == pd.Timestamp("2019-01-06T12:00:00")
    sample_bar = _progress_text("Download", 5242880, 10485760, "bytes")
    assert "50.0%" in sample_bar and "5.0/10.0 MiB" in sample_bar
    assert "total size unknown" in _progress_text("Download", 524288, None, "bytes")
    with TemporaryDirectory() as tmp:
        output = Path(tmp) / "test_candidates.csv"
        record = dict(start=t0, end=t0 + pd.Timedelta(minutes=10),
                      core_start=t0 + pd.Timedelta(minutes=5),
                      core_end=t0 + pd.Timedelta(minutes=6),
                      category="screened_full", peak_angle_deg=89.9,
                      screening_method="synthetic")
        first = _save_scan_checkpoint(
            output, t0, t0 + pd.Timedelta(days=2), t0.date(),
            [record], 86400
        )
        assert output.is_file() and len(first) == 1
        saved_day, saved_rows, valid = _load_scan_checkpoint(
            output, t0, t0 + pd.Timedelta(days=2)
        )
        assert saved_day == t0 and valid == 86400
        assert len(saved_rows) == 1
        assert saved_rows[0]["core_start"] == record["core_start"]
        assert len(pd.read_csv(output)) == 1
        try:
            _load_scan_checkpoint(output, t0, t0 + pd.Timedelta(days=3))
        except ValueError:
            pass
        else:
            raise AssertionError("Checkpoint range validation failed")

    # A mid-stream disconnect must keep the partial and request only the rest.
    from types import SimpleNamespace
    from unittest.mock import patch
    with TemporaryDirectory() as tmp:
        target = Path(tmp) / "test_download.cdf"
        seen_headers = []

        class FakeTransfer:
            def __init__(self, status, headers, data, fail=False):
                self.status_code = status
                self.headers = headers
                self.data = data
                self.fail = fail

            def __enter__(self):
                return self

            def __exit__(self, *args):
                return False

            def raise_for_status(self):
                pass

            def iter_content(self, size):
                yield self.data
                if self.fail:
                    raise requests.exceptions.ChunkedEncodingError(
                        "Synthetic incomplete NASA download")

        def interrupted_get(url, **kwargs):
            seen_headers.append(kwargs['headers'])
            if len(seen_headers) == 1:
                return FakeTransfer(
                    200, {'Content-Length': '2048'}, b'A'*1024, fail=True)
            return FakeTransfer(
                206, {'Content-Length': '1024',
                      'Content-Range': 'bytes 1024-2047/2048'},
                b'B'*1024)

        with patch("requests.get", side_effect=interrupted_get), \
             patch("cdflib.CDF") as cdf_open:
            cdf_open.return_value.cdf_info.return_value = SimpleNamespace(
                zVariables=[FULL_VECTOR, FULL_EPOCH], rVariables=[])
            output = _download_full_cdf(
                "https://example.org/synthetic.cdf", target, attempts=2,
                pause=lambda seconds: None)
        assert output.read_bytes() == b'A'*1024 + b'B'*1024
        assert 'Range' not in seen_headers[0]
        assert seen_headers[1]['Range'] == 'bytes=1024-'
        assert not target.with_suffix('.cdf.part').exists()

    print("PASS: screening, progress, incremental save/resume, and range parsing")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Read and pre-screen Parker Solar Probe full-cadence magnetic data"
    )
    parser.add_argument("--start", default="2018-11-06",
                        help="UTC start date/time, e.g. 2018-11-06T03:30:00")
    parser.add_argument("--end", default=None,
                        help="UTC end date/time (exclusive)")
    parser.add_argument("--duration", default=None,
                        help="Duration from --start: 10m, 2h, 5d, 3w, 2mo")
    parser.add_argument("--days", type=int, default=1,
                        help="Default number of days if --end/--duration not provided")
    parser.add_argument("--resolution", choices=tuple(PRODUCTS), default="1min",
                        help="Resolution for direct reading without --scan")
    parser.add_argument("--scan", "--screen-full", dest="scan", action="store_true",
                        help="Automatically pre-screen the entire selected time range")
    parser.add_argument("--output", type=Path, default=None,
                        help="Optional location for the candidate-window CSV")
    parser.add_argument("--resume", action="store_true",
                        help="Continue an interrupted scan from its saved checkpoint")
    parser.add_argument("--self-test", action="store_true",
                        help="Test preliminary screening without NASA downloads")
    args = parser.parse_args()
    if args.self_test:
        _screen_self_test()
    elif args.scan:
        if args.end is not None and args.duration is not None:
            parser.error("Choose either --end or --duration")
        start = _as_utc_naive(args.start)
        if args.duration is not None:
            end = _parse_duration(start, args.duration)
        elif args.end is not None:
            end = _as_utc_naive(args.end)
        else:
            if args.days < 1:
                parser.error("--days must be positive")
            end = start + pd.Timedelta(days=args.days)
        def report(day, candidates, valid, expected):
            print(f"[{day}] {valid}/{expected} valid one-second bins; "
                  f"{candidates} preliminary windows (before merging)", flush=True)
        output = args.output or (
            Path("data") /
            f"magnetic_scan_{start:%Y%m%dT%H%M%S}_{end:%Y%m%dT%H%M%S}.csv"
        )
        if args.resume:
            try:
                resume_after, previous_rows, prior_valid = _load_scan_checkpoint(
                    output, start, end
                )
            except (FileNotFoundError, ValueError) as exc:
                parser.error(str(exc))
            print(f"[RESUME] Continuing after {resume_after.date()} "
                  f"with {len(previous_rows)} original candidates", flush=True)
        else:
            resume_after, previous_rows, prior_valid = None, [], 0
        def save_day(day, rows, valid):
            _save_scan_checkpoint(output, start, end, day, rows, valid)
        selection = scan_time_range(
            start, end, on_day=report, on_checkpoint=save_day,
            resume_after=resume_after, initial_records=previous_rows,
            initial_valid_seconds=prior_valid
        )
        # Retain the completed scan checkpoint as a verified coverage
        # manifest. --keep-strong-only in magnetic_data_analysis.py needs
        # its requested UTC bounds and valid-second coverage before
        # authorizing irreversible raw-cache cleanup.
        _atomic_text(output, selection.to_csv(index=False))
        print(f"COMPLETE: {len(selection)} investigation windows; saved {output}")
        print(f"[COVERAGE MANIFEST] {_checkpoint_file(output)}", flush=True)
    else:
        if args.end is not None or args.duration is not None:
            parser.error("--end/--duration require --scan for time-range screening")
        if args.days < 1:
            parser.error("--days must be positive")
        for day in days(args.start, args.days):
            frame = get_day(day, resolution=args.resolution)
            print(f"{day} [{args.resolution}]: {frame.Bmag.notna().sum()} "
                  f"valid samples of {len(frame)}")
