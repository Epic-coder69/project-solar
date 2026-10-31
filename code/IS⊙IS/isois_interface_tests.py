#!/usr/bin/env python3
"""Offline regression tests for Project Solar's unified IS⊙IS CLI routing.

Tests command dispatch and preservation of legacy aliases without reading,
downloading, rewriting or deleting any research data.
Run: python code/IS⊙IS/isois_interface_tests.py
"""
from __future__ import annotations

import importlib.util
import io
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
IS_DIR = Path(__file__).resolve().parent


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise AssertionError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


reader = load(IS_DIR / "isois_reader.py", "isois_reader_test_target")
analyzer = load(IS_DIR / "isois_analyzer.py", "isois_analyzer_test_target")
project = load(ROOT / "project.py", "project_test_target")


class UnifiedInterfaces(unittest.TestCase):
    def test_reader_epihi_dispatch_without_implicit_download(self):
        with patch.object(reader.subprocess, "call", return_value=0) as call:
            rc = reader.main(["epihi", "--events-csv", "catalog.csv",
                              "--resume", "--products", "second-rates"])
        self.assertEqual(rc, 0)
        argv = call.call_args.args[0]
        self.assertEqual(argv[0], sys.executable)
        self.assertTrue(argv[1].endswith("particle_highres_reader.py"))
        self.assertEqual(argv[2:], ["--events-csv", "catalog.csv", "--resume",
                                     "--products", "second-rates"])
        self.assertNotIn("--download", argv)

    def test_summary_dispatch(self):
        with patch.object(reader.subprocess, "call", return_value=0) as call:
            self.assertEqual(reader.main(["summary", "--offline",
                                          "--file", "x.cdf"]), 0)
        self.assertTrue(call.call_args.args[0][1].endswith("particle_data_reader.py"))
        self.assertEqual(call.call_args.args[0][2:],
                         ["--offline", "--file", "x.cdf"])

    def test_unknown_mode_does_not_run_backend(self):
        with patch.object(reader.subprocess, "call") as call:
            with redirect_stdout(io.StringIO()):
                self.assertEqual(reader.main(["unknown"]), 2)
        call.assert_not_called()

    def test_reader_help_never_invokes_backend(self):
        with patch.object(reader.subprocess, "call") as call:
            with redirect_stdout(io.StringIO()):
                self.assertEqual(reader.main(["--help"]), 0)
        call.assert_not_called()

    def test_legacy_highres_option_dispatch(self):
        with patch.object(reader.subprocess, "call", return_value=0) as call:
            self.assertEqual(reader.main(["--highres", "--plan-only"]), 0)
        self.assertEqual(call.call_args.args[0][2:], ["--plan-only"])

    def test_coverage_preserves_explicit_catalog(self):
        with patch.object(analyzer.subprocess, "call", return_value=0) as call:
            self.assertEqual(analyzer.main(
                ["coverage", "--events-csv", "my.csv",
                 "--before-seconds", "10"]), 0)
        args = call.call_args.args[0]
        self.assertTrue(args[1].endswith("particle_multires_coverage.py"))
        self.assertEqual(args.count("--events-csv"), 1)
        self.assertEqual(args[args.index("--events-csv") + 1], "my.csv")
        self.assertNotIn("--download", args)

    def test_analyzer_modes_forward_to_distinct_existing_backends(self):
        for mode, name in (
            ("fine-coverage", "particle_coverage_audit.py"),
            ("survey", "particle_data_analysis.py"),
        ):
            with self.subTest(mode=mode):
                with patch.object(analyzer.subprocess, "call",
                                  return_value=0) as call:
                    self.assertEqual(
                        analyzer.main([mode, "--help"]), 0)
                self.assertTrue(call.call_args.args[0][1].endswith(name))

    def test_analyzer_uses_local_default_catalog_only_if_present(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            file = root / analyzer.EVENTS
            file.parent.mkdir(parents=True)
            file.write_text("event_id\n")
            with patch.object(analyzer, "ROOT", root):
                args = analyzer._prepared_args("coverage", [])
                self.assertEqual(args, ["--events-csv", str(analyzer.EVENTS)])
                self.assertEqual(
                    analyzer._prepared_args("coverage",
                                            ["--events-csv=custom.csv"]),
                    ["--events-csv=custom.csv"])
                self.assertEqual(
                    analyzer._prepared_args("coverage", ["--self-test"]),
                    ["--self-test"])

    def test_project_legacy_aliases(self):
        aliases = {
            "particles": ("isois_reader.py", "summary"),
            "highres": ("isois_reader.py", "epihi"),
            "particle-analysis": ("isois_analyzer.py", "survey"),
            "fine-coverage": ("isois_analyzer.py", "fine-coverage"),
            "coverage": ("isois_analyzer.py", "coverage"),
        }
        for command, (backend, mode) in aliases.items():
            with self.subTest(command=command):
                with patch.object(project.subprocess, "call",
                                  return_value=0) as call:
                    self.assertEqual(
                        project.invoke(command, ["--help"]), 0)
                argv = call.call_args.args[0]
                self.assertTrue(argv[1].endswith(backend))
                self.assertEqual(argv[2:], [mode, "--help"])
                self.assertEqual(call.call_args.kwargs["cwd"], str(project.ROOT))

    def test_new_entrypoints_and_data_paths_unchanged(self):
        self.assertTrue(project.SCRIPTS["isois-reader"].endswith("isois_reader.py"))
        self.assertTrue(project.SCRIPTS["isois-analyzer"].endswith("isois_analyzer.py"))
        self.assertEqual(
            project.REPORT,
            Path("data/isois/highres/multires_coverage_summary.json"))
        self.assertEqual(
            project.EVENTS,
            Path("data/full_analysis_nov_dec_2018_candidates_strong_events.csv"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
