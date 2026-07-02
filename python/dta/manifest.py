#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 11:31:10 2026

@author: jesse
"""
from __future__ import annotations

from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Any, Iterator, Literal
from dta.utils import validate_path, infer_leaflet


Leaflet = Literal["upper", "lower"]
FileType = Literal["raw counts", "avg counts"]


@dataclass(frozen=True, slots=True)
class DTAFileRecord:
    """
    A single file that participates in a DTA analysis.

    Arguments
    ---------
    path : Path | str
        The path to the file. Stored as Path after validation.

    Attributes
    ----------
    leaflet : Leaflet
        "upper" or "lower". Inferred from filename.
    filetype : FileType
        "raw counts" or "avg counts". Inferred from filename.
    user_str1, user_str2, user_str3 : str
        User-defined fields, each containing a string that can be used to slice
        records in the DTAProjectManifest.
    user_int1: int
        User-defined field containing an int that can be used to slice records
        in the DTAProjectManifest.
    """

    path: Path | str
    leaflet: Leaflet = field(init=False)
    filetype: FileType = field(init=False)

    user_str1: str | None = None
    user_str2: str | None = None
    user_str3: str | None = None
    user_int1: int | None = None
    user_int2: int | None = None

    def __post_init__(self) -> None:
        """Validate path and infer leaflet and filetype."""
        object.__setattr__(self, "path", validate_path(self.path))
        if self.path.suffix != ".dat":
            raise ValueError(f"path should end in .dat, not {self.path.suffix}.")
        object.__setattr__(self, "leaflet", infer_leaflet(self.path))
        if ".avg." in self.path.name:
            object.__setattr__(self, "filetype", "avg counts")
        else:
            object.__setattr__(self, "filetype", "raw counts")

    def check(self) -> None:
        """Raise an error if the file does not exist."""
        if not self.path.is_file():
            raise FileNotFoundError(f"No file found at {self.path}.")


@dataclass(slots=True)
class DTAProjectManifest:
    """
    A curated collection of DTA file records.

    Arguments
    ---------
    records : list[DTAFileRecord], Optional
        User may supply records on creation, and/or can use the ``add`` method later.

    Attributes
    ----------
    records : list[DTAFileRecord]
        The list of DTAFileRecords contained in this Manifest.
    record_fields : set[str]
        A set of the attributes available in a DTAFileRecord.

    """

    records: list[DTAFileRecord] = field(default_factory=list)

    def __post_init__(self) -> None:
        """Run user-supplied records through self.add for validation."""
        records = self.records
        self.records = []
        for record in records:
            self.add(record)

    def __len__(self) -> int:
        """Return length of self.records."""
        return len(self.records)

    def __iter__(self) -> Iterator[DTAFileRecord]:
        """Iterate through self.records."""
        return iter(self.records)

    @property
    def record_fields(self) -> set[str]:
        """List all the attributes inside of a DTAFileRecord."""
        return {f.name for f in fields(DTAFileRecord)}

    def add(self, record: DTAFileRecord) -> None:
        """Add a DTAFileRecord to self.records."""
        if not isinstance(record, DTAFileRecord):
            raise TypeError(f"Expected DTAFileRecord, not {type(record)}.")
        if any(existing.path == record.path for existing in self.records):
            raise ValueError(f"{record.path} already exists in the manifest.")
        self.records.append(record)

    def select(self, **filters: Any) -> DTAProjectManifest:
        """
        Return a new manifest containing records that satisfy the specified
        filters.
    
        Each keyword argument specifies a field of ``DTAFileRecord`` and one or
        more acceptable values for that field.
    
        Filters applied to different fields are combined using a logical **AND**.
        Therefore, a record must satisfy **every** specified filter to be included
        in the returned manifest.
    
        If the value associated with a filter is a ``list``, ``tuple``, or ``set``,
        the elements of that collection are combined using a logical **OR**. A
        record satisfies that filter if its field value matches **any** value in
        the collection.
    
        Parameters
        ----------
        **filters
            Keyword arguments specifying one or more ``DTAFileRecord`` fields and
            their desired values. Filter names must correspond to fields of
            ``DTAFileRecord``.
    
        Returns
        -------
        DTAProjectManifest
            A new manifest containing the records that satisfy the specified
            filters. The original manifest is not modified.
    
        Raises
        ------
        ValueError
            If no filters are supplied.
        ValueError
            If one or more filter names are not fields of ``DTAFileRecord``.
    
        Examples
        --------
        Select all records from the upper leaflet.
    
        >>> upper = manifest.select(leaflet="upper")
    
        Select records whose ``user_str1`` is either ``"PC_CL"`` or
        ``"PC_CHOL"``.
    
        >>> selected = manifest.select(
        ...     user_str1={"PC_CL", "PC_CHOL"},
        ... )
    
        Select records from the upper leaflet **and** whose ``user_str1`` is
        either ``"PC_CL"`` **or** ``"PC_CHOL"``. This corresponds to the logical
        expression::
    
            leaflet == "upper" AND
            (user_str1 == "PC_CL" OR user_str1 == "PC_CHOL")
    
        >>> selected = manifest.select(
        ...     leaflet="upper",
        ...     user_str1={"PC_CL", "PC_CHOL"},
        ... )
        """
        if not filters:
            raise ValueError("At least one filter must be supplied.")
        self._validate_record_fields(filters.keys())
        new_record_list = []
        for record in self.records:
            if self._record_matches(record, filters):
                new_record_list.append(record)
        return DTAProjectManifest(new_record_list)

    def _validate_record_fields(self, field_names: Any) -> None:
        """Prevent mispelled DTAFileRecord attributes from causing problems."""
        unknown = set(field_names) - self.record_fields
        if unknown:
            raise ValueError(
                f"Unknown DTAFileRecord field(s): {sorted(unknown)}. "
                f"Valid fields are: {sorted(self.record_fields)}."
            )

    @staticmethod
    def _record_matches(record: DTAFileRecord, filters: dict[str, Any]) -> bool:
        """
        Determine whether a record satisfies a set of filters.

        Each key in ``filters`` must be the name of a ``DTAFileRecord`` attribute.
        Filters applied to different attributes are combined using a logical
        **AND**. If the value associated with a filter is a ``list``, ``tuple``,
        or ``set``, the elements of that collection are combined using a logical
        **OR**.

        Parameters
        ----------
        record : DTAFileRecord
            The record to test.
        filters : dict[str, Any]
            Dictionary mapping ``DTAFileRecord`` attribute names to one or more
            acceptable values.

        Returns
        -------
        bool
            ``True`` if the record satisfies every filter and ``False`` otherwise.

        Notes
        -----
        This method assumes that the filter names have already been validated by
        :meth:`_validate_record_fields`.
        """
        for field_name, expected in filters.items():
            actual = getattr(record, field_name)

            if isinstance(expected, (set, list, tuple)):
                if actual not in expected:
                    return False
            elif actual != expected:
                return False

        return True
