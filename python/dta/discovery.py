#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 19:47:29 2026

@author: jesse
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

from dta.manifest import DTAFileRecord


def discover_files(
    folders: Path | str | Iterable[Path | str],
    *,
    recursive: bool = True,
) -> list[DTAFileRecord]:
    """
    Discover DTA ``.dat`` files in one or more folders.

    This function scans the supplied folder or folders, creates a
    ``DTAFileRecord`` for each discovered ``.dat`` file, and returns the
    successfully constructed records.

    Parameters
    ----------
    folders : Path | str | Iterable[Path | str]
        One folder, or an iterable of folders, to search.
    recursive : bool, default=True
        If True, search each folder recursively. If False, search only the
        immediate contents of each folder.

    Returns
    -------
    list[DTAFileRecord]
        Records created from discovered ``.dat`` files.

    Raises
    ------
    FileNotFoundError
        If any supplied folder does not exist.
    NotADirectoryError
        If any supplied path is not a directory.
    """

    folder_list = _convert_dir_to_path(folders)

    records: list[DTAFileRecord] = []

    for folder in folder_list:
        _validate_folder(folder)

        pattern = "**/*.dat" if recursive else "*.dat"

        for path in sorted(folder.glob(pattern)):
            try:
                record = DTAFileRecord(path)
                records.append(record)
            except Exception:
                continue
    return records


def _convert_dir_to_path(
    folders: Path | str | Iterable[Path | str],
) -> list[Path]:
    """Convert one or more directories to a list of ``Path`` objects."""
    if isinstance(folders, (str, Path)):
        return [Path(folders)]
    return [Path(folder) for folder in folders]


def _validate_folder(folder: Path) -> None:
    """Raise an error if ``folder`` is not an existing directory."""
    if not folder.exists():
        raise FileNotFoundError(f"No folder found at {folder}.")
    if not folder.is_dir():
        raise NotADirectoryError(f"{folder} is not a directory.")
