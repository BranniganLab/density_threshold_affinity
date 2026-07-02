"""Unit tests for DTA file discovery utilities.

These tests exercise the public ``discover_files`` function and the helper
functions that normalize and validate input folders.
"""

from pathlib import Path

import pytest

from dta.discovery import (
    _convert_dir_to_path,
    _validate_folder,
    discover_files,
)
from dta.manifest import DTAFileRecord


def touch(path: Path) -> Path:
    """Create an empty file and return its path."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("")
    return path


def test_convert_dir_to_path_accepts_single_string_path(tmp_path):
    """Test that a single string path is converted to a one-item Path list."""
    result = _convert_dir_to_path(str(tmp_path))

    assert result == [tmp_path]


def test_convert_dir_to_path_accepts_single_path(tmp_path):
    """Test that a single Path is converted to a one-item Path list."""
    result = _convert_dir_to_path(tmp_path)

    assert result == [tmp_path]


def test_convert_dir_to_path_accepts_iterable_of_paths_and_strings(tmp_path):
    """Test that an iterable of path-like inputs is converted to Path objects."""
    folder1 = tmp_path / "folder1"
    folder2 = tmp_path / "folder2"

    result = _convert_dir_to_path([folder1, str(folder2)])

    assert result == [folder1, folder2]


def test_validate_folder_accepts_existing_directory(tmp_path):
    """Test that _validate_folder accepts an existing directory."""
    _validate_folder(tmp_path)


def test_validate_folder_raises_for_missing_directory(tmp_path):
    """Test that _validate_folder raises FileNotFoundError for a missing path."""
    missing = tmp_path / "missing"

    with pytest.raises(FileNotFoundError, match="No folder found"):
        _validate_folder(missing)


def test_validate_folder_raises_for_file_instead_of_directory(tmp_path):
    """Test that _validate_folder raises NotADirectoryError for a file path."""
    file_path = touch(tmp_path / "not_a_directory.dat")

    with pytest.raises(NotADirectoryError, match="is not a directory"):
        _validate_folder(file_path)


def test_discover_files_finds_valid_dat_files_in_single_folder(tmp_path):
    """Test that discover_files returns DTAFileRecord objects for valid files."""
    upper = touch(tmp_path / "CL.upp.dat")
    lower = touch(tmp_path / "CL.low.dat")

    records = discover_files(tmp_path, recursive=False)

    assert len(records) == 2
    assert all(isinstance(record, DTAFileRecord) for record in records)
    assert {record.path for record in records} == {upper, lower}
    assert {record.leaflet for record in records} == {"upper", "lower"}


def test_discover_files_accepts_multiple_folders(tmp_path):
    """Test that discover_files can search more than one input folder."""
    folder1 = tmp_path / "dataset1"
    folder2 = tmp_path / "dataset2"

    file1 = touch(folder1 / "CL.upp.dat")
    file2 = touch(folder2 / "CHOL.low.dat")

    records = discover_files([folder1, folder2], recursive=False)

    assert {record.path for record in records} == {file1, file2}


def test_discover_files_recursive_true_finds_nested_files(tmp_path):
    """Test that recursive discovery finds files in nested subdirectories."""
    top_level = touch(tmp_path / "CL.upp.dat")
    nested = touch(tmp_path / "rep1" / "CL.low.dat")

    records = discover_files(tmp_path, recursive=True)

    assert {record.path for record in records} == {top_level, nested}


def test_discover_files_recursive_false_ignores_nested_files(tmp_path):
    """Test that non-recursive discovery ignores files below the top folder."""
    top_level = touch(tmp_path / "CL.upp.dat")
    touch(tmp_path / "rep1" / "CL.low.dat")

    records = discover_files(tmp_path, recursive=False)

    assert len(records) == 1
    assert records[0].path == top_level


def test_discover_files_ignores_non_dat_files(tmp_path):
    """Test that only files ending in .dat are considered."""
    valid = touch(tmp_path / "CL.upp.dat")
    touch(tmp_path / "notes.txt")
    touch(tmp_path / "CL.upp.csv")

    records = discover_files(tmp_path, recursive=False)

    assert len(records) == 1
    assert records[0].path == valid


def test_discover_files_skips_invalid_dat_files(tmp_path):
    """Test that invalid .dat files are skipped by the current implementation.

    DTAFileRecord raises when it cannot infer the leaflet from the filename.
    discover_files currently catches that exception and continues.
    """
    valid = touch(tmp_path / "CL.upp.dat")
    touch(tmp_path / "invalid_name.dat")

    records = discover_files(tmp_path, recursive=False)

    assert len(records) == 1
    assert records[0].path == valid


def test_discover_files_returns_records_in_sorted_path_order(tmp_path):
    """Test that records are returned in deterministic sorted path order."""
    b_file = touch(tmp_path / "B.low.dat")
    a_file = touch(tmp_path / "A.upp.dat")

    records = discover_files(tmp_path, recursive=False)

    assert [record.path for record in records] == [a_file, b_file]


def test_discover_files_raises_for_missing_input_folder(tmp_path):
    """Test that discover_files validates input folders before scanning."""
    missing = tmp_path / "missing"

    with pytest.raises(FileNotFoundError, match="No folder found"):
        discover_files(missing)


def test_discover_files_raises_for_file_input(tmp_path):
    """Test that discover_files rejects file paths as input folders."""
    file_path = touch(tmp_path / "CL.upp.dat")

    with pytest.raises(NotADirectoryError, match="is not a directory"):
        discover_files(file_path)
