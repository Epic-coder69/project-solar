#!/usr/bin/env python3
"""Unified IS⊙IS analysis entry point for Project Solar.

Separates scientific/audit steps from the CDF ingestion interface.
All output filenames and CSV/JSON schemas are currently those of the tested
existing analysis backends; this router does not rewrite their measurements.

Modes:
  coverage      — preferred-channel multiresolution *timing* audit
  fine-coverage — exact native 1–10 s FIELDS event/flank timing audit
  survey        — exploratory hourly count-rate survey characterization

Quality-flag interpretation, calibrated passband harmonization and Poisson
significance of particle acceleration are NOT implemented by these audits.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EVENTS = Path("data/full_analysis_nov_dec_2018_candidates_strong_events.csv")
SURVEY = Path("data/isois/fields_event_summary_rates.csv")
BACKENDS = {
    "coverage": "particle_multires_coverage.py",
    "fine-coverage": "particle_coverage_audit.py",
    "survey": "particle_data_analysis.py",
}


def help_text() -> None:
    print("""Project Solar | Unified IS⊙IS analyzer

Usage:
  python project.py isois-analyzer coverage [coverage options]
  python project.py isois-analyzer fine-coverage [fine coverage options]
  python project.py isois-analyzer survey --rates-csv FILE --events-csv FILE
  python project.py isois-analyzer --self-test

  python project.py isois-analyzer coverage --help
  python project.py isois-analyzer fine-coverage --help
  python project.py isois-analyzer survey --help

When running a coverage mode without --events-csv, the existing 599-candidate
catalog in data/ is used if present. No new data is downloaded.
All results retain their existing locations, filenames and schemas.

Coverage reports show NOMINAL timestamp overlap, not instrument-quality
verification, accurate cross-product energy equivalence, or particle-
acceleration significance. Hourly observations never count as fine timing.
""")


def _has_flag(args: list[str], flag: str) -> bool:
    return flag in args or any(arg.startswith(flag + "=") for arg in args)


def _prepared_args(mode: str, args: list[str]) -> list[str]:
    """Choose optional local defaults without changing explicit user inputs."""
    prepared = list(args)
    if mode in {"coverage", "fine-coverage"} and not (
        _has_flag(prepared, "--events-csv") or
        "--self-test" in prepared or
        "-h" in prepared or
        "--help" in prepared
    ):
        if (ROOT / EVENTS).is_file():
            prepared.extend(["--events-csv", str(EVENTS)])
    if mode == "coverage" and not (
        _has_flag(prepared, "--summary-csv") or
        "--self-test" in prepared or
        "-h" in prepared or
        "--help" in prepared
    ):
        # This file is optional: multires analyzer will use saved EPI-Hi
        # hourly partitions even if the separate summary CSV is unavailable.
        if (ROOT / SURVEY).is_file():
            prepared.extend(["--summary-csv", str(SURVEY)])
    return prepared


def _call(mode: str, args: list[str]) -> int:
    backend = HERE / BACKENDS[mode]
    if not backend.is_file():
        print(f"[ERROR] Missing IS⊙IS analysis backend: {backend}",
              file=sys.stderr)
        return 2
    prepared = _prepared_args(mode, args)
    return subprocess.call([sys.executable, str(backend), *prepared],
                           cwd=str(ROOT))


def self_test() -> int:
    for mode in ("fine-coverage", "coverage", "survey"):
        print(f"[ANALYZER TEST] {mode}", flush=True)
        rc = _call(mode, ["--self-test"])
        if rc:
            return rc
    # Verify default injection is disabled for explicit flags and self-tests.
    assert _prepared_args("coverage", ["--self-test"]) == ["--self-test"]
    assert _prepared_args("coverage", ["--events-csv=custom.csv"])[
        :1] == ["--events-csv=custom.csv"]
    assert _prepared_args("survey", ["--rates-csv", "x"]) == [
        "--rates-csv", "x"]
    print("PASS: unified IS⊙IS analyzer dispatches existing scientific audits")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in {"help", "--help", "-h"}:
        help_text()
        return 0
    if args == ["--self-test"]:
        return self_test()
    mode, *remaining = args
    if mode not in BACKENDS:
        print(f"[ERROR] Unknown IS⊙IS analyzer mode: {mode}", file=sys.stderr)
        help_text()
        return 2
    return _call(mode, remaining)


if __name__ == "__main__":
    raise SystemExit(main())
