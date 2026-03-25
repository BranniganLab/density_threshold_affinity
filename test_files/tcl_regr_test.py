#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 15:59:07 2026

@author: js2746
"""

from pathlib import Path
import difflib

def assert_directories_equal(dir_a: Path, dir_b: Path, *, binary=False):
    files_a = {p.relative_to(dir_a) for p in dir_a.rglob("*") if p.is_file()}
    files_b = {p.relative_to(dir_b) for p in dir_b.rglob("*") if p.is_file()}

    # Check file set equality first
    missing = files_a - files_b
    extra = files_b - files_a

    if missing or extra:
        msg = []
        if missing:
            msg.append(f"Missing in {dir_b}:\n" + "\n".join(map(str, sorted(missing))))
        if extra:
            msg.append(f"Extra in {dir_b}:\n" + "\n".join(map(str, sorted(extra))))
        raise AssertionError("\n\n".join(msg))

    # Compare file contents
    for rel_path in sorted(files_a):
        file_a = dir_a / rel_path
        file_b = dir_b / rel_path

        if binary:
            assert file_a.read_bytes() == file_b.read_bytes(), f"Binary files differ: {rel_path}"
        else:
            a_lines = file_a.read_text().splitlines()
            b_lines = file_b.read_text().splitlines()

            if a_lines != b_lines:
                diff = "\n".join(
                    difflib.unified_diff(
                        a_lines, b_lines,
                        fromfile=str(file_a),
                        tofile=str(file_b),
                        lineterm=""
                    )
                )
                raise AssertionError(f"Files differ: {rel_path}\n{diff}")


def test_output_directory():
    assert_directories_equal(
        Path("./MD_files/rep1/ref_vals"),
        Path("./MD_files/rep1/test_vals"),
    )
