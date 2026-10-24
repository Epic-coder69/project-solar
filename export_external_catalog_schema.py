#!/usr/bin/env python3
"""Export a small, reviewable schema report for Martinovic et al. switchback HDF5.

Does not commit or duplicate the underlying third-party measurements.
Run after: python external_switchback_validation.py --inspect
"""
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

SOURCE = Path("data/external/ApJ_SB_Database.h5")
DEST = Path("results/external_catalog_schema.json")


def compact(value, limit=4):
    """JSON-safe sample, bounded in size; preserve structured field names."""
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="replace")[:180]
    if isinstance(value, np.void) and value.dtype.names:
        return {k: compact(value[k], limit) for k in value.dtype.names}
    if isinstance(value, np.ndarray):
        arr = value.reshape(-1)
        return [compact(x, limit) for x in arr[:limit]]
    if isinstance(value, np.generic):
        return compact(value.item(), limit)
    if isinstance(value, (tuple, list)):
        return [compact(x, limit) for x in value[:limit]]
    if isinstance(value, (str, float, int, bool)) or value is None:
        return str(value)[:180] if isinstance(value, str) else value
    return str(value)[:180]


def main():
    if not SOURCE.exists():
        raise SystemExit(f"Missing {SOURCE}; run external_switchback_validation.py --inspect first.")
    report = {
        "source": "https://doi.org/10.5281/zenodo.4579821",
        "filename": SOURCE.name,
        "md5": hashlib.md5(SOURCE.read_bytes()).hexdigest(),
        "published_md5": "5046a6964a2f17383c243892202ad5c5",
        "objects": [],
        "pandas_tables": [],
    }
    if report["md5"] != report["published_md5"]:
        raise SystemExit("STOP: checksum differs from published HDF5.")
    with h5py.File(SOURCE, "r") as h5:
        def visit(name, item):
            entry = {"path": "/" + name, "type": "dataset" if isinstance(item, h5py.Dataset) else "group"}
            entry["attributes"] = {
                str(k): compact(v) for k, v in list(item.attrs.items())[:15]
            }
            if isinstance(item, h5py.Dataset):
                entry["shape"] = list(item.shape)
                entry["dtype"] = str(item.dtype)
                entry["fields"] = list(item.dtype.names or [])
                try:
                    if item.ndim == 0:
                        sample = item[()]
                    elif item.shape[0] == 0:
                        sample = []
                    else:
                        sample = item[: min(item.shape[0], 2)]
                    entry["sample"] = compact(sample)
                except Exception as exc:
                    entry["sample_error"] = str(exc)[:160]
            report["objects"].append(entry)
        h5.visititems(visit)

    # Some HDF5 files are written by pandas/PyTables. Show those keys separately.
    try:
        import pandas as pd
        with pd.HDFStore(SOURCE, mode="r") as store:
            for key in store.keys():
                entry = {"key": key}
                try:
                    df = store.select(key, start=0, stop=2)
                    entry["columns"] = [str(c) for c in df.columns]
                    entry["sample"] = [compact(v) for v in df.to_dict(orient="records")]
                except Exception as exc:
                    entry["preview_error"] = str(exc)[:160]
                report["pandas_tables"].append(entry)
    except Exception as exc:
        report["pandas_inspection_note"] = str(exc)[:250]

    DEST.parent.mkdir(parents=True, exist_ok=True)
    DEST.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("Saved", DEST)
    print("HDF5 objects:", len(report["objects"]))
    print("Pandas tables:", len(report["pandas_tables"]))
    print("Only a small schema report is committed; raw HDF5 remains in data/.")


if __name__ == "__main__":
    main()
