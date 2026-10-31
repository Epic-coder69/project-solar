#!/usr/bin/env python3
"""Unified IS⊙IS data-ingestion interface for Project Solar.

One public reader with two deliberately distinct NASA products:
    epihi   — native 1s/10s/60s/300s/3600s EPI-Hi CDF -> NPZ
    summary — original L2 summary count-rate survey -> CSV

These subcommands delegate to the existing tested decoders. No measurements
are resampled, renamed as identical detector channels, combined or deleted.
Legacy reader entry points remain available for reproducibility.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BACKENDS = {
    "epihi": "particle_highres_reader.py",
    "summary": "particle_data_reader.py",
}


def help_text() -> None:
    print("""Project Solar | Unified IS⊙IS reader

Usage:
  python project.py isois-reader epihi   [EPI-Hi reader options]
  python project.py isois-reader summary [summary reader options]
  python project.py isois-reader list-products
  python project.py isois-reader --self-test

Examples:
  python project.py isois-reader epihi --events-csv \
data/full_analysis_nov_dec_2018_candidates_strong_events.csv \
--products let1-rates10 het-rates10 second-rates --download --resume

  python project.py isois-reader summary --help
  python project.py isois-reader epihi --help

EPI-Hi products preserve source-specific integration times, energy bins,
counts, calibrated flux and measurement quality flags where available.
Summary records remain a separate hourly count-rate survey: do NOT treat
them as event-resolved or as equivalent to EPI-Hi differential flux.
No NASA downloads occur without an explicit --download flag.
All existing NPZ/CSV paths, manifest formats and cache directories remain
unchanged. Older reader scripts and project.py aliases still work.
""")


def _call(backend: str, args: list[str]) -> int:
    """Isolate legacy argument parsers and preserve their exit statuses."""
    path = HERE / BACKENDS[backend]
    if not path.is_file():
        print(f"[ERROR] Missing IS⊙IS decoder: {path}", file=sys.stderr)
        return 2
    return subprocess.call([sys.executable, str(path), *args], cwd=str(HERE.parents[1]))


def list_products() -> int:
    from particle_highres_reader import CHANNEL_GUIDE, PRODUCTS, DEFAULT_PRODUCTS

    print("NASA IS⊙IS L2 SUMMARY: separate hourly survey rates")
    print("NASA EPI-Hi native-cadence products:")
    for name in PRODUCTS:
        default = " [default]" if name in DEFAULT_PRODUCTS else ""
        print(f"  {name}{default}: {CHANNEL_GUIDE[name]}")
    print("Use --products with the epihi subcommand to select products.")
    return 0


def self_test() -> int:
    """Run actual existing offline decoder regressions, not mock success."""
    for mode in ("summary", "epihi"):
        print(f"[READER TEST] {mode}", flush=True)
        rc = _call(mode, ["--self-test"])
        if rc:
            return rc
    print("PASS: unified IS⊙IS reader routes both existing decoders")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in {"help", "-h", "--help"}:
        help_text()
        return 0
    if args == ["--self-test"]:
        return self_test()
    if args[0] in {"list-products", "--list-products"}:
        if len(args) != 1:
            print("[ERROR] list-products takes no arguments", file=sys.stderr)
            return 2
        return list_products()
    # Compatibility with old --highres option in direct reader invocations.
    if args[0] == "--highres":
        args = ["epihi", *args[1:]]
    if args[0] == "highres":
        args[0] = "epihi"
    mode, *remaining = args
    if mode not in BACKENDS:
        print(f"[ERROR] Unknown IS⊙IS reader mode: {mode}", file=sys.stderr)
        help_text()
        return 2
    return _call(mode, remaining)


if __name__ == "__main__":
    raise SystemExit(main())
