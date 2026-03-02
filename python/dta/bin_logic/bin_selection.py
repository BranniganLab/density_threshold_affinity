#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Selection model for discrete 2D bin indices.

This module defines the `BinSelection` class, a minimal state container
for managing a set of selected bins represented as `(r_index, theta_index)`
integer tuples. It provides controlled mutation and retrieval of the
selection state and is intended to serve as the core selection model
within the DTA GUI architecture.

The class currently supports full replacement, retrieval, and clearing
of the selection. Undo/redo hooks are reserved for future implementation.
"""

from __future__ import annotations
from collections.abc import Iterable
from typing import Union
from dta.bin_logic import BinAddress, BinAddressLike, as_bin_address


class BinSelection:
    """Class that contains the set of selected bins."""

    def __init__(self) -> None:
        """Initialize an empty selection."""
        self._bins: set[BinAddress] = set()

    def get_bins(self) -> set[BinAddress]:
        """Return a copy of the current selection."""
        return set(self._bins)

    def set_bins(self, bins: Iterable[BinAddressLike]) -> None:
        """Replace/set the current selection."""
        self._bins = {as_bin_address(b) for b in bins}

    def clear(self) -> None:
        """Clear the selection."""
        self._bins.clear()

    def undo(self):
        """Future undo function."""

    def redo(self):
        """Future redo function."""
