#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 14:15:52 2026.

@author: js2746
"""


class BinSelection:
    """Class that contains the set of selected bins."""

    def __init__(self):
        """Initialize an empty selection."""
        self._bins = set()

    def get_bins(self):
        """
        Return the current selection.

        Returns
        -------
        set[tuple[int, int]]
        """
        return set(self._bins)

    def set_bins(self, bins):
        """
        Replace/set the current selection.

        Parameters
        ----------
        bins : iterable of (int, int)
            New selection.
        """
        self._bins = set(bins)

    def clear(self):
        """Clear the selection."""
        self._bins.clear()

    def undo(self):
        """Future undo function."""

    def redo(self):
        """Future redo function."""
