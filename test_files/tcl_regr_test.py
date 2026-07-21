#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Regression test comparing the contents of two directory trees.

Created on Wed Mar 25 15:59:07 2026

@author: js2746
"""

from difflib import SequenceMatcher
from itertools import zip_longest
from pathlib import Path
import fnmatch

import pytest


def is_ignored(rel_path: Path, patterns: list[str]) -> bool:
    """
    Determine whether a relative file path should be ignored.

    A path is considered ignored if it matches any pattern in ``patterns``.
    Matching is performed against both the full relative path and the filename.

    Parameters
    ----------
    rel_path : pathlib.Path
        Path relative to the root directory being compared.
    patterns : list[str]
        Filename, relative-path, or glob patterns to ignore.

    Returns
    -------
    bool
        ``True`` if the path matches an ignore pattern; otherwise ``False``.
    """
    rel_str = str(rel_path)
    name = rel_path.name

    return any(
        fnmatch.fnmatch(rel_str, pattern)
        or fnmatch.fnmatch(name, pattern)
        for pattern in patterns
    )


def _format_line_difference(
    expected_line: str | None,
    actual_line: str | None,
    line_number: int,
) -> list[str]:
    """Format one differing pair of lines with field-level differences."""
    if expected_line is None:
        return [
            f"  line {line_number}: unexpected line in actual file",
            f"    actual:   {actual_line}",
        ]

    if actual_line is None:
        return [
            f"  line {line_number}: line missing from actual file",
            f"    expected: {expected_line}",
        ]

    expected_fields = expected_line.split()
    actual_fields = actual_line.split()

    output = [f"  line {line_number}:"]

    if len(expected_fields) == len(actual_fields):
        differing_fields = [
            (index, expected, actual)
            for index, (expected, actual) in enumerate(
                zip(expected_fields, actual_fields),
                start=1,
            )
            if expected != actual
        ]

        if differing_fields:
            for index, expected, actual in differing_fields:
                output.append(
                    f"    field {index}: expected {expected!r}; actual {actual!r}"
                )
            return output

    # Field-level alignment is ambiguous when the number of fields differs,
    # or when only whitespace differs. Show the complete lines instead.
    output.extend(
        [
            f"    expected: {expected_line}",
            f"    actual:   {actual_line}",
        ]
    )
    return output


def _format_text_file_difference(
    expected_file: Path,
    actual_file: Path,
    rel_path: Path,
) -> str | None:
    """Return a detailed report for a text mismatch, or ``None`` if equal."""
    expected_lines = expected_file.read_text().splitlines()
    actual_lines = actual_file.read_text().splitlines()

    if expected_lines == actual_lines:
        return None

    report = [f"TEXT MISMATCH: {rel_path}"]
    mismatch_count = 0
    matcher = SequenceMatcher(a=expected_lines, b=actual_lines, autojunk=False)

    for tag, expected_start, expected_end, actual_start, actual_end in matcher.get_opcodes():
        if tag == "equal":
            continue

        expected_region = expected_lines[expected_start:expected_end]
        actual_region = actual_lines[actual_start:actual_end]

        for offset, (expected_line, actual_line) in enumerate(
            zip_longest(expected_region, actual_region)
        ):
            mismatch_count += 1
            if expected_line is None:
                line_number = actual_start + offset + 1
            else:
                line_number = expected_start + offset + 1

            report.extend(
                _format_line_difference(expected_line, actual_line, line_number)
            )

    report.insert(
        1,
        f"  {mismatch_count} differing line(s); "
        f"expected {len(expected_lines)} line(s), actual {len(actual_lines)} line(s)",
    )
    return "\n".join(report)


def assert_directories_equal(
    dir_a: Path,
    dir_b: Path,
    *,
    binary: bool = False,
    ignore: list[str] | None = None,
) -> None:
    """
    Assert that two directories have identical structure and file contents.

    All structural and content mismatches are collected before the assertion
    fails. For text files, differences are reported by line and whitespace-
    separated field so changed values are easy to identify.

    Parameters
    ----------
    dir_a, dir_b : pathlib.Path
        Reference and test-output directories, respectively.
    binary : bool, optional
        Compare files as raw bytes instead of text when ``True``.
    ignore : list[str], optional
        Filename, relative-path, or glob patterns to ignore.

    Raises
    ------
    pytest.fail.Exception
        If files are missing, extra, or differ in content.
    """
    ignore = ignore or []

    files_a = {
        path.relative_to(dir_a)
        for path in dir_a.rglob("*")
        if path.is_file()
        and not is_ignored(path.relative_to(dir_a), ignore)
    }
    files_b = {
        path.relative_to(dir_b)
        for path in dir_b.rglob("*")
        if path.is_file()
        and not is_ignored(path.relative_to(dir_b), ignore)
    }

    missing = sorted(files_a - files_b)
    extra = sorted(files_b - files_a)
    common = sorted(files_a & files_b)
    mismatch_reports: list[str] = []

    if missing:
        mismatch_reports.append(
            f"MISSING FROM ACTUAL DIRECTORY ({dir_b}):\n  "
            + "\n  ".join(map(str, missing))
        )

    if extra:
        mismatch_reports.append(
            f"EXTRA IN ACTUAL DIRECTORY ({dir_b}):\n  "
            + "\n  ".join(map(str, extra))
        )

    for rel_path in common:
        expected_file = dir_a / rel_path
        actual_file = dir_b / rel_path

        if binary:
            if expected_file.read_bytes() != actual_file.read_bytes():
                mismatch_reports.append(f"BINARY MISMATCH: {rel_path}")
            continue

        try:
            report = _format_text_file_difference(
                expected_file,
                actual_file,
                rel_path,
            )
        except UnicodeDecodeError as error:
            mismatch_reports.append(
                f"TEXT DECODE ERROR: {rel_path}\n"
                f"  Could not decode as text: {error}"
            )
            continue

        if report is not None:
            mismatch_reports.append(report)

    if mismatch_reports:
        mismatch_count = len(mismatch_reports)
        summary = (
            "Directory comparison failed with "
            f"{mismatch_count} mismatch group(s):"
        )
        report = summary + "\n\n" + "\n\n".join(mismatch_reports)

        # Keep the detailed diagnostics out of the exception message. Pytest
        # otherwise repeats a multiline exception in both the failure traceback
        # and the short test summary.
        print(report, flush=True)
        pytest.fail(
            f"Directory comparison found {mismatch_count} mismatch group(s).",
            pytrace=False,
        )


def test_output_directory():
    """Test that reference and generated output directories match."""
    assert_directories_equal(
        Path("./MD_files/rep1/ref_vals"),
        Path("./MD_files/rep1/test_vals"),
        ignore=[".placeholder"],
    )
