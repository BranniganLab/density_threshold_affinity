"""Tests for DTAFileRecord and DTAProjectManifest.

These tests cover the manifest data model, including file-record construction,
automatic metadata inference, duplicate-path validation, and generic record
selection using user-defined fields.
"""

from pathlib import Path

import pytest

from dta.manifest import DTAFileRecord, DTAProjectManifest


def make_dat_file(tmp_path: Path, name: str) -> Path:
    """Create a minimal .dat file and return its path."""
    path = tmp_path / name
    path.write_text("header\n", encoding="utf-8")
    return path


@pytest.fixture
def records(tmp_path: Path) -> list[DTAFileRecord]:
    """Create a representative set of records for manifest tests."""
    return [
        DTAFileRecord(
            make_dat_file(tmp_path, "CL.upp.dat"),
            user_str1="PC_CL",
            user_str2="CL",
            user_str3="run_A",
            user_int1=1,
            user_int2=10,
        ),
        DTAFileRecord(
            make_dat_file(tmp_path, "CL.low.dat"),
            user_str1="PC_CL",
            user_str2="CL",
            user_str3="run_A",
            user_int1=2,
            user_int2=10,
        ),
        DTAFileRecord(
            make_dat_file(tmp_path, "CHOL.upp.avg.dat"),
            user_str1="PC_CHOL",
            user_str2="CHOL",
            user_str3="run_B",
            user_int1=1,
            user_int2=20,
        ),
    ]


def test_file_record_stores_path_as_path_object(tmp_path: Path) -> None:
    """DTAFileRecord converts string paths to pathlib.Path objects."""
    path = make_dat_file(tmp_path, "CL.upp.dat")

    record = DTAFileRecord(str(path))

    assert isinstance(record.path, Path)
    assert record.path == path


def test_file_record_infers_upper_leaflet_from_filename(tmp_path: Path) -> None:
    """DTAFileRecord infers upper leaflet from filenames containing '.upp.'."""
    record = DTAFileRecord(make_dat_file(tmp_path, "CL.upp.dat"))

    assert record.leaflet == "upper"


def test_file_record_infers_lower_leaflet_from_filename(tmp_path: Path) -> None:
    """DTAFileRecord infers lower leaflet from filenames containing '.low.'."""
    record = DTAFileRecord(make_dat_file(tmp_path, "CL.low.dat"))

    assert record.leaflet == "lower"


def test_file_record_infers_avg_counts_filetype(tmp_path: Path) -> None:
    """DTAFileRecord labels filenames containing '.avg.' as average counts."""
    record = DTAFileRecord(make_dat_file(tmp_path, "CL.upp.avg.dat"))

    assert record.filetype == "avg counts"


def test_file_record_infers_raw_counts_filetype(tmp_path: Path) -> None:
    """DTAFileRecord labels non-average .dat files as raw counts."""
    record = DTAFileRecord(make_dat_file(tmp_path, "CL.upp.dat"))

    assert record.filetype == "raw counts"


def test_file_record_rejects_non_dat_file(tmp_path: Path) -> None:
    """DTAFileRecord raises ValueError when the path does not end in .dat."""
    path = tmp_path / "CL.upp.txt"
    path.write_text("header\n", encoding="utf-8")

    with pytest.raises(ValueError, match="path should end in .dat"):
        DTAFileRecord(path)


def test_file_record_rejects_missing_file(tmp_path: Path) -> None:
    """DTAFileRecord raises FileNotFoundError when the file does not exist."""
    with pytest.raises(FileNotFoundError):
        DTAFileRecord(tmp_path / "missing.upp.dat")


def test_file_record_check_raises_if_file_removed(tmp_path: Path) -> None:
    """check raises FileNotFoundError if a previously valid file is removed."""
    path = make_dat_file(tmp_path, "CL.upp.dat")
    record = DTAFileRecord(path)
    path.unlink()

    with pytest.raises(FileNotFoundError):
        record.check()


def test_manifest_len_and_iter(records: list[DTAFileRecord]) -> None:
    """DTAProjectManifest supports len() and iteration over records."""
    manifest = DTAProjectManifest(records)

    assert len(manifest) == 3
    assert list(manifest) == records


def test_manifest_add_appends_record(records: list[DTAFileRecord]) -> None:
    """add appends a DTAFileRecord to the manifest."""
    manifest = DTAProjectManifest()

    manifest.add(records[0])

    assert len(manifest) == 1
    assert manifest.records[0] is records[0]


def test_manifest_add_rejects_non_record() -> None:
    """add raises TypeError for objects that are not DTAFileRecord instances."""
    manifest = DTAProjectManifest()

    with pytest.raises(TypeError, match="Expected DTAFileRecord"):
        manifest.add("not a record")  # type: ignore[arg-type]


def test_manifest_add_rejects_duplicate_path(records: list[DTAFileRecord]) -> None:
    """add raises ValueError when a record path is already present."""
    manifest = DTAProjectManifest([records[0]])

    with pytest.raises(ValueError, match="already exists in the manifest"):
        manifest.add(records[0])


def test_manifest_constructor_rejects_duplicate_paths(
    records: list[DTAFileRecord],
) -> None:
    """Constructor validates user-supplied records through add()."""
    with pytest.raises(ValueError, match="already exists in the manifest"):
        DTAProjectManifest([records[0], records[0]])


def test_record_fields_includes_record_dataclass_fields(
    records: list[DTAFileRecord],
) -> None:
    """record_fields exposes all dataclass fields available for selection."""
    manifest = DTAProjectManifest(records)

    assert {
        "path",
        "leaflet",
        "filetype",
        "user_str1",
        "user_str2",
        "user_str3",
        "user_int1",
        "user_int2",
    }.issubset(manifest.record_fields)


def test_select_with_single_filter(records: list[DTAFileRecord]) -> None:
    """select returns records matching one field value."""
    manifest = DTAProjectManifest(records)

    selected = manifest.select(leaflet="upper")

    assert len(selected) == 2
    assert all(record.leaflet == "upper" for record in selected)


def test_select_combines_different_fields_with_logical_and(
    records: list[DTAFileRecord],
) -> None:
    """select combines filters on different fields using logical AND."""
    manifest = DTAProjectManifest(records)

    selected = manifest.select(user_str1="PC_CL", leaflet="upper")

    assert len(selected) == 1
    assert selected.records[0].user_str1 == "PC_CL"
    assert selected.records[0].leaflet == "upper"


def test_select_combines_collection_values_with_logical_or(
    records: list[DTAFileRecord],
) -> None:
    """select treats collection values within one filter as logical OR."""
    manifest = DTAProjectManifest(records)

    selected = manifest.select(user_str2={"CL", "CHOL"})

    assert len(selected) == 3


def test_select_combines_and_between_fields_and_or_within_field(
    records: list[DTAFileRecord],
) -> None:
    """select applies AND between filters and OR within collection filters."""
    manifest = DTAProjectManifest(records)

    selected = manifest.select(
        leaflet="upper",
        user_str1={"PC_CL", "PC_CHOL"},
        user_int1=1,
    )

    assert len(selected) == 2
    assert all(record.leaflet == "upper" for record in selected)
    assert all(record.user_int1 == 1 for record in selected)


def test_select_returns_new_manifest_without_modifying_original(
    records: list[DTAFileRecord],
) -> None:
    """select returns a new manifest and leaves the original unchanged."""
    manifest = DTAProjectManifest(records)

    selected = manifest.select(user_str1="PC_CL")

    assert selected is not manifest
    assert len(selected) == 2
    assert len(manifest) == 3


def test_select_with_no_filters_raises_value_error():
    """
    Test that calling select() without any filters raises a ValueError.
    """
    manifest = DTAProjectManifest()

    with pytest.raises(
        ValueError,
        match="At least one filter must be supplied",
    ):
        manifest.select()


def test_select_rejects_unknown_filter_name(records: list[DTAFileRecord]) -> None:
    """select raises ValueError when a filter is not a DTAFileRecord field."""
    manifest = DTAProjectManifest(records)

    with pytest.raises(ValueError, match="Unknown DTAFileRecord field"):
        manifest.select(misspelled_field="PC_CL")


def test_record_matches_returns_true_when_all_filters_match(
    records: list[DTAFileRecord],
) -> None:
    """_record_matches returns True when a record satisfies all filters."""
    record = records[0]

    assert DTAProjectManifest._record_matches(
        record,
        {"leaflet": "upper", "user_str1": "PC_CL"},
    )


def test_record_matches_returns_false_when_any_filter_fails(
    records: list[DTAFileRecord],
) -> None:
    """_record_matches returns False when any required filter fails."""
    record = records[0]

    assert not DTAProjectManifest._record_matches(
        record,
        {"leaflet": "lower", "user_str1": "PC_CL"},
    )


def test_record_matches_supports_allowed_value_collections(
    records: list[DTAFileRecord],
) -> None:
    """_record_matches treats list, tuple, and set values as allowed values."""
    record = records[0]

    assert DTAProjectManifest._record_matches(record, {"user_str1": ["PC_CL"]})
    assert DTAProjectManifest._record_matches(record, {"user_str1": ("PC_CL",)})
    assert DTAProjectManifest._record_matches(record, {"user_str1": {"PC_CL"}})
