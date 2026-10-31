#!/usr/bin/env python3
"""Project Solar entry point. Run from any directory with this script on disk.

Examples:
  python project.py status
  python project.py test
  python project.py coverage
  python project.py highres --events-csv data/...
  python project.py sweap --self-test

IS⊙IS now exposes one public reader and one public analyzer with the
original tested instrument-specific backends retained underneath. Legacy
command aliases and all stored NASA data paths are preserved.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = {
    "fields-reader": "code/FIELDS/magnetic_data_reader.py",
    "fields-analysis": "code/FIELDS/magnetic_data_analysis.py",
    "isois-reader": "code/IS⊙IS/isois_reader.py",
    "isois-analyzer": "code/IS⊙IS/isois_analyzer.py",
    # Legacy command names remain valid via the same unified interfaces.
    "particles": "code/IS⊙IS/isois_reader.py",
    "highres": "code/IS⊙IS/isois_reader.py",
    "particle-analysis": "code/IS⊙IS/isois_analyzer.py",
    "fine-coverage": "code/IS⊙IS/isois_analyzer.py",
    "coverage": "code/IS⊙IS/isois_analyzer.py",
    "sweap": "code/SWEAP/sweap_data_reader.py",
}
SELF_TESTS = [
    "fields-reader",
    "fields-analysis",
    "isois-reader",
    "isois-analyzer",
    "sweap",
]
EVENTS = Path("data/full_analysis_nov_dec_2018_candidates_strong_events.csv")
SURVEY = Path("data/isois/fields_event_summary_rates.csv")
REPORT = Path("data/isois/highres/multires_coverage_summary.json")


def help_text() -> None:
    print("""Project Solar — stable command interface

  python project.py status                    Check local data and report paths
  python project.py test                      Run all eight offline instrument self-tests
  python project.py isois-reader --help       Unified IS⊙IS CDF reader
  python project.py isois-analyzer --help     Unified particle coverage and survey analysis
  python project.py isois-reader list-products
  python project.py isois-reader epihi --help EPI-Hi products/cadences
  python project.py isois-analyzer coverage   Recalculate multires timing coverage
  python project.py highres --help            Legacy EPI-Hi command (still works)
  python project.py particles --help          Legacy hourly survey command
  python project.py coverage                  Legacy multires coverage command
  python project.py fine-coverage --help
  python project.py particle-analysis --help
  python project.py fields-reader --help
  python project.py fields-analysis --help
  python project.py sweap --help              SWEAP CDF inspection

Subcommands forward their options to the existing tested backends. No NASA
data downloads occur unless --download is explicitly requested.
Git tracks the source code, not the local data/ directory.
""")


def status() -> int:
    print("Project Solar directory:", ROOT)
    bad = False
    for label, path in (
        ("FIELDS catalog", EVENTS),
        ("IS⊙IS summary CSV", SURVEY),
        ("Multires coverage report", REPORT),
    ):
        exists = (ROOT / path).is_file()
        print(f"  {'FOUND' if exists else 'MISSING'} {label}: {path}")
    for label, relative in SCRIPTS.items():
        if not (ROOT / relative).is_file():
            print(f"  MISSING code for {label}: {relative}")
            bad = True
    for old in ("code/M_FIELD", "code/IS_O_IS"):
        if (ROOT / old).exists():
            print(f"  LEGACY FOLDER {old} (local leftover; not tracked by Git)")
    print("NOTE: MISSING local research files are not necessarily errors;")
    print("      data/ is Git-ignored and is not downloaded by git pull.")
    return 1 if bad else 0


def invoke(command: str, args: list[str]) -> int:
    source = ROOT / SCRIPTS[command]
    if not source.is_file():
        print(f"Missing script: {source}", file=sys.stderr)
        return 2
    legacy_modes = {
        "particles": "summary",
        "highres": "epihi",
        "particle-analysis": "survey",
        "fine-coverage": "fine-coverage",
        "coverage": "coverage",
    }
    if command in legacy_modes:
        args = [legacy_modes[command], *args]
    return subprocess.call(
        [sys.executable, str(source), *args], cwd=str(ROOT)
    )


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in {"-h", "--help", "help"}:
        help_text()
        return 0
    command, *args = argv
    if command == "status":
        if args:
            print("status takes no arguments", file=sys.stderr)
            return 2
        return status()
    if command == "test":
        if args:
            print("test takes no arguments", file=sys.stderr)
            return 2
        for name in SELF_TESTS:
            print(f"\n[TEST] {name}", flush=True)
            rc = invoke(name, ["--self-test"])
            if rc:
                print(f"[FAILED] {name} exit={rc}", file=sys.stderr)
                return rc
        # Exercise the new command routing without NASA access or writes.
        routing = ROOT / "code/IS⊙IS/isois_interface_tests.py"
        if not routing.is_file():
            print(f"[FAILED] Missing IS⊙IS routing tests: {routing}",
                  file=sys.stderr)
            return 2
        print("\n[TEST] IS⊙IS unified-interface routing", flush=True)
        rc = subprocess.call([sys.executable, str(routing)], cwd=str(ROOT))
        if rc:
            print(f"[FAILED] IS⊙IS routing tests exit={rc}", file=sys.stderr)
            return rc
        print("\nPASS: all 8 instrument self-tests and unified routing tests")
        return 0
    if command not in SCRIPTS:
        print(f"Unknown command: {command}", file=sys.stderr)
        help_text()
        return 2
    return invoke(command, args)


if __name__ == "__main__":
    raise SystemExit(main())
