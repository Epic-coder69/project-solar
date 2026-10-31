#!/usr/bin/env python3
"""Project Solar entry point. Run from any directory with this script on disk.

Examples:
  python project.py status
  python project.py test
  python project.py coverage
  python project.py highres --events-csv data/...
  python project.py sweap --self-test

Scripts retain their original locations; this CLI avoids Unicode-path and
instrument-folder renames in user commands. Large NASA data stay local.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = {
    "fields-reader": "code/FIELDS/magnetic_data_reader.py",
    "fields-analysis": "code/FIELDS/magnetic_data_analysis.py",
    "particles": "code/IS⊙IS/particle_data_reader.py",
    "highres": "code/IS⊙IS/particle_data_reader.py",
    "particle-analysis": "code/IS⊙IS/particle_data_analysis.py",
    "fine-coverage": "code/IS⊙IS/particle_coverage_audit.py",
    "coverage": "code/IS⊙IS/particle_multires_coverage.py",
    "sweap": "code/SWEAP/sweap_data_reader.py",
}
SELF_TESTS = [
    "fields-reader",
    "fields-analysis",
    "particles",
    "particle-analysis",
    "highres",
    "fine-coverage",
    "coverage",
    "sweap",
]
EVENTS = Path("data/full_analysis_nov_dec_2018_candidates_strong_events.csv")
SURVEY = Path("data/isois/fields_event_summary_rates.csv")
REPORT = Path("data/isois/highres/multires_coverage_summary.json")


def help_text() -> None:
    print("""Project Solar — local Python entry point

  python project.py status               Check code and local data/report paths
  python project.py cleanup              Preview duplicate folders/logs safely
  python project.py cleanup --apply      Remove confirmed duplicates and logs
  python project.py test                 Run all eight offline synthetic self-tests
  python project.py coverage             Recalculate 599-event multires timing audit
  python project.py fine-coverage --help Inspect 10-second flank coverage options
  python project.py highres --help       Native EPI-Hi acquisition and extraction
  python project.py particles --help     Legacy hourly summary survey reader
  python project.py particle-analysis --help
  python project.py fields-reader --help
  python project.py fields-analysis --help
  python project.py sweap --help         SWEAP CDF inspection

Arguments after a command are forwarded unchanged to the original script.
No command downloads NASA data unless you explicitly use its --download
option. Git tracks source code, not the local data/ directory.
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
            print(f"  LEGACY FOLDER {old} — preview: python project.py cleanup")
    if (ROOT / "logs").exists():
        print("  LOCAL LOGS logs/ — preview: python project.py cleanup")
    print("NOTE: MISSING local research files are not necessarily errors;")
    print("      data/ is Git-ignored and is not downloaded by git pull.")
    return 1 if bad else 0


def invoke(command: str, args: list[str]) -> int:
    source = ROOT / SCRIPTS[command]
    if not source.is_file():
        print(f"Missing script: {source}", file=sys.stderr)
        return 2
    if command == "highres":
        args = ["--highres", *args]
    if command == "coverage" and not args:
        if not (ROOT / EVENTS).is_file():
            print(f"Missing local event catalog: {EVENTS}", file=sys.stderr)
            return 2
        args = ["--events-csv", str(EVENTS)]
        if (ROOT / SURVEY).is_file():
            args += ["--summary-csv", str(SURVEY)]
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
    if command == "cleanup":
        return subprocess.call(
            [sys.executable, str(ROOT / "cleanup_local.py"), *args],
            cwd=str(ROOT)
        )
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
        print("\nPASS: all 8 offline instrument self-tests")
        return 0
    if command not in SCRIPTS:
        print(f"Unknown command: {command}", file=sys.stderr)
        help_text()
        return 2
    return invoke(command, args)


if __name__ == "__main__":
    raise SystemExit(main())
