#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Persistent selection state for polar bin grids.

This module defines the domain-level selection model used to store and
manipulate the set of selected bins independently of any GUI framework.

The selection model:
- represents committed selection state only (no transient previews),
- supports replace/add/remove operations at the bin level,
- is consumed by both GUI controllers and analysis code.

It contains no rendering logic and no assumptions about user input
mechanisms. All interaction semantics are handled by GUI layers, and all
visualization is handled by renderers.
"""


class BinSelection:
    """
    Mutable domain model representing a set of selected bins.

    This class contains no rendering logic and no interaction logic.
    All selection mutations are centralized here to enable future
    undo/redo support.
    """

    def __init__(self):
        """Initialize an empty selection."""
        self._bins = set()

    def snapshot(self):
        """
        Capture the current selection state.

        Returns
        -------
        frozenset
            Immutable snapshot of selected bins.
        """
        return frozenset(self._bins)

    def update(self, bins):
        """
        Replace the current selection.

        Parameters
        ----------
        bins : iterable of (int, int)
            New selection.
        """
        self._bins = set(bins)

    def clear(self):
        """Clear the selection."""
        self._bins.clear()

    def get_bins(self):
        """
        Return the current selection.

        Returns
        -------
        set[tuple[int, int]]
        """
        return set(self._bins)
