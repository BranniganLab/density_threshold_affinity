#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 15:59:07 2026

@author: js2746
"""

from pathlib import Path
import difflib
import fnmatch


def is_ignored(rel_path: Path, patterns: list[str]) -> bool:
    """
    Determine whether a relative file path should be ignored.

    A path is considered ignored if it matches any pattern in `patterns`.
    Matching is performed against both:
      - the full relative path (e.g., "subdir/file.txt")
      - the filename only (e.g., "file.txt")

    Supported pattern types (via fnmatch):
      - Exact filenames: "file.txt"
      - Relative paths: "subdir/file.txt"
      - Glob patterns: "*.log", "cache/*", "*/tmp/*"

    Parameters
    ----------
    rel_path : Path
        Path relative to the root directory being compared.
    patterns : list[str]
        List of ignore patterns.

    Returns
    -------
    bool
        True if the path matches any ignore pattern, False otherwise.
    """
    rel_str = str(rel_path)
    name = rel_path.name

    for pattern in patterns:
        if (
            fnmatch.fnmatch(rel_str, pattern)
            or fnmatch.fnmatch(name, pattern)
        ):
            return True

    return False


def assert_directories_equal(
    dir_a: Path,
    dir_b: Path,
    *,
    binary: bool = False,
    ignore: list[str] | None = None,
):
    """
    Assert that two directories are identical in structure and file contents.

    Parameters
    ----------
    dir_a, dir_b : Path
        Directories to compare.
    binary : bool, optional
        If True, compare files as raw bytes.
    ignore : list[str], optional
        Patterns to ignore. Supports:
        - filenames: "log.txt"
        - relative paths: "subdir/file.txt"
        - glob patterns: "*.log", "*/tmp/*"
    """
    ignore = ignore or []

    files_a = {
        p.relative_to(dir_a)
        for p in dir_a.rglob("*")
        if p.is_file() and not is_ignored(p.relative_to(dir_a))
    }

    files_b = {
        p.relative_to(dir_b)
        for p in dir_b.rglob("*")
        if p.is_file() and not is_ignored(p.relative_to(dir_b))
    }

    # --- structure check ---
    missing = files_a - files_b
    extra = files_b - files_a

    if missing or extra:
        msg = []
        if missing:
            msg.append(
                f"Missing in {dir_b}:\n" + "\n".join(map(str, sorted(missing)))
            )
        if extra:
            msg.append(
                f"Extra in {dir_b}:\n" + "\n".join(map(str, sorted(extra)))
            )
        raise AssertionError("\n\n".join(msg))

    # --- content check ---
    for rel_path in sorted(files_a):
        file_a = dir_a / rel_path
        file_b = dir_b / rel_path

        if binary:
            if file_a.read_bytes() != file_b.read_bytes():
                raise AssertionError(f"Binary files differ: {rel_path}")
        else:
            a_lines = file_a.read_text().splitlines()
            b_lines = file_b.read_text().splitlines()

            if a_lines != b_lines:
                diff = "\n".join(
                    difflib.unified_diff(
                        a_lines,
                        b_lines,
                        fromfile=str(file_a),
                        tofile=str(file_b),
                        lineterm="",
                    )
                )
                raise AssertionError(
                    f"Files differ: {rel_path}\n{diff}"
                )


def test_output_directory():
    """Test if contents of reference val and test val directories match."""
    assert_directories_equal(
        Path("./MD_files/rep1/ref_vals"),
        Path("./MD_files/rep1/test_vals"),
        ignore=[".placeholder"],
    )
