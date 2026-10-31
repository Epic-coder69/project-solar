#!/usr/bin/env python3
"""Safely remove obsolete local instrument-folder duplicates and disposable logs.

Default: preview only. Actual deletion requires --apply.
Never edits or deletes the canonical instrument code or data/ science files.
"""
from __future__ import annotations

import argparse
import filecmp
import os
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ALIASES = (
    ("code/M_FIELD", "code/FIELDS"),
    ("code/IS_O_IS", "code/IS⊙IS"),
)
CACHE_DIRS = {"__pycache__", ".pytest_cache"}
TRASH_FILES = {".DS_Store", ".gitkeep"}


@dataclass(frozen=True)
class Action:
    kind: str       # file or directory tree
    path: Path
    why: str


def _safe_file(path: Path) -> bool:
    return path.is_file() and not path.is_symlink()


def _safe_directory(path: Path) -> bool:
    return path.is_dir() and not path.is_symlink()


def plan(root: Path) -> tuple[list[Action], list[str], list[Path]]:
    """Allow-list based: only obsolete aliases, root logs, root *.log files."""
    root = root.resolve(strict=True)
    planned: list[Action] = []
    warnings: list[str] = []
    old_directories: list[Path] = []

    for old_name, current_name in ALIASES:
        old = root / old_name
        current = root / current_name
        if old.is_symlink():
            warnings.append(f"Skipping symlinked legacy directory: {old_name}")
            continue
        if not old.exists():
            continue
        if not _safe_directory(old):
            warnings.append(f"Skipping unexpected legacy path: {old_name}")
            continue
        if not _safe_directory(current):
            warnings.append(f"Skipping {old_name}: current {current_name} missing or unsafe")
            continue

        old_directories.append(old)
        for folder, dirs, files in os.walk(old, topdown=True, followlinks=False):
            here = Path(folder)
            remaining_dirs = []
            for d in dirs:
                path = here / d
                if path.is_symlink():
                    warnings.append(f"Preserving symlink: {path.relative_to(root)}")
                elif d in CACHE_DIRS:
                    planned.append(Action("tree", path, "rebuildable Python cache"))
                else:
                    remaining_dirs.append(d)
            dirs[:] = remaining_dirs
            for name in files:
                path = here / name
                rel = path.relative_to(old)
                replacement = current / rel
                if path.is_symlink():
                    warnings.append(f"Preserving symlink: {path.relative_to(root)}")
                    continue
                if not _safe_file(path):
                    warnings.append(f"Preserving unexpected entry: {path.relative_to(root)}")
                    continue
                if name in TRASH_FILES or path.suffix in {".pyc", ".pyo"}:
                    planned.append(Action("file", path, "disposable legacy metadata/cache"))
                elif _safe_file(replacement) and filecmp.cmp(
                    path, replacement, shallow=False
                ):
                    planned.append(Action(
                        "file", path,
                        f"byte-identical copy of {replacement.relative_to(root)}"
                    ))
                else:
                    warnings.append(
                        f"PRESERVED unique/different file: {path.relative_to(root)}"
                    )

    # The user requested removal of generated log output. No data/ files
    # (including any diagnostic files there) are considered by this utility.
    logs = root / "logs"
    if logs.is_symlink():
        warnings.append("Preserving symlinked logs/ path")
    elif _safe_directory(logs):
        planned.append(Action("tree", logs, "obsolete local logs folder"))
    elif _safe_file(logs):
        planned.append(Action("file", logs, "obsolete local logs file"))
    for file in root.glob("*.log"):
        if _safe_file(file):
            planned.append(Action("file", file, "generated root log"))
        elif file.is_symlink():
            warnings.append(f"Preserving symlinked log: {file.name}")

    return planned, warnings, old_directories


def cleanup(root: Path, *, apply: bool) -> tuple[int, int]:
    root = root.resolve(strict=True)
    actions, warnings, old_directories = plan(root)
    print(f"[CLEANUP] {'APPLY' if apply else 'PREVIEW ONLY'}")
    print("Allowed: identical files in code/M_FIELD or code/IS_O_IS;")
    print("         old Python caches, root logs/ and root *.log files.")
    print("Protected: data/, canonical code/FIELDS and code/IS⊙IS,")
    print("           different/unique files, original CDFs and NPZs.")
    for action in actions[:50]:
        print(f"  {'REMOVE' if apply else 'WOULD REMOVE'} "
              f"{action.path.relative_to(root)} [{action.why}]")
    if len(actions) > 50:
        print(f"  ... plus {len(actions) - 50} more matching items")
    for message in warnings:
        print("  WARNING:", message)
    if not apply:
        print(f"[DRY RUN] {len(actions)} possible items; nothing deleted.")
        print("To perform these specific deletions: python project.py cleanup --apply")
        return len(actions), len(warnings)

    removed = 0
    for action in actions:
        try:
            if action.path.is_symlink():
                print(f"  SKIP changed-to-symlink: {action.path.relative_to(root)}")
                continue
            if action.kind == "tree" and _safe_directory(action.path):
                shutil.rmtree(action.path)
                removed += 1
            elif action.kind == "file" and _safe_file(action.path):
                # Recheck byte-identity at deletion time if it was the reason.
                if action.why.startswith("byte-identical copy of "):
                    current = root / action.why.removeprefix("byte-identical copy of ")
                    if not _safe_file(current) or not filecmp.cmp(
                        action.path, current, shallow=False
                    ):
                        print(f"  SKIP changed file: {action.path.relative_to(root)}")
                        continue
                action.path.unlink()
                removed += 1
        except OSError as exc:
            print(f"  WARNING: could not remove {action.path.relative_to(root)}: {exc}")

    for old in old_directories:
        # Remove only genuinely empty legacy directories. If any unique
        # work remains, it is preserved with the folder.
        for folder, dirs, files in os.walk(old, topdown=False, followlinks=False):
            location = Path(folder)
            if location.is_symlink():
                continue
            try:
                location.rmdir()
                print(f"  REMOVED empty legacy directory: {location.relative_to(root)}")
            except OSError:
                pass
        if old.exists():
            print(f"  PRESERVED legacy folder with remaining files: {old.relative_to(root)}")
    print(f"[COMPLETE] Removed {removed} matching items; "
          f"{len(warnings)} warnings/preserved entries.")
    return removed, len(warnings)


def self_test() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        for name in ("code/FIELDS", "code/M_FIELD", "code/IS⊙IS",
                     "code/IS_O_IS", "data/isois/highres", "logs"):
            (root / name).mkdir(parents=True)
        (root / "code/FIELDS/magnetic_data_reader.py").write_text("same")
        (root / "code/M_FIELD/magnetic_data_reader.py").write_text("same")
        (root / "code/IS⊙IS/particle_data_reader.py").write_text("current")
        (root / "code/IS_O_IS/particle_data_reader.py").write_text("different")
        (root / "code/IS_O_IS/local_experiment.py").write_text("unique")
        cache = root / "code/IS_O_IS/__pycache__"
        cache.mkdir()
        (cache / "reader.pyc").write_bytes(b"cache")
        (root / "logs/execution.txt").write_text("disposable")
        (root / "session.log").write_text("disposable")
        science = root / "data/isois/highres/measurements.npz"
        science.write_bytes(b"important science")

        n, _ = cleanup(root, apply=False)
        assert n >= 4
        assert (root / "logs").exists()
        assert (root / "code/M_FIELD").exists()
        n, _ = cleanup(root, apply=True)
        assert n >= 4
        assert not (root / "logs").exists()
        assert not (root / "session.log").exists()
        assert not (root / "code/M_FIELD").exists()
        assert not cache.exists()
        assert science.read_bytes() == b"important science"
        assert (root / "code/IS⊙IS/particle_data_reader.py").read_text() == "current"
        assert (root / "code/IS_O_IS/local_experiment.py").exists()
        assert (root / "code/IS_O_IS/particle_data_reader.py").exists()
        count, _ = cleanup(root, apply=True)
        assert count == 0
    print("PASS: preview, identical duplicates, logs removal, unique-file "
          "preservation and science data protection")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--apply", action="store_true",
                       help="actually remove safe items (default: preview)")
    group.add_argument("--self-test", action="store_true",
                       help="run disposable-fixture regression test")
    args = parser.parse_args(argv)
    if args.self_test:
        self_test()
    else:
        cleanup(ROOT, apply=args.apply)


if __name__ == "__main__":
    main()
