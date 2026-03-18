#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for BinSelection state management and input normalization."""

from dta.bin_logic import BinAddress, BinSelection


def test_bin_selection_initially_empty():
    """Verify that a new BinSelection starts with no selected bins."""
    sel = BinSelection()
    assert sel.get_bins() == set()


def test_bin_selection_set_bins_replaces_selection_and_normalizes_duplicates():
    """Verify that set_bins replaces prior contents and collapses duplicate inputs."""
    sel = BinSelection()

    # Verify tuple inputs are normalized and duplicates collapse to unique BinAddress values.
    sel.set_bins([(0, 1), (2, 3), (0, 1)])
    assert sel.get_bins() == {BinAddress(0, 1), BinAddress(2, 3)}

    # Verify a later call replaces, rather than merges with, the prior selection.
    sel.set_bins([(9, 9)])
    assert sel.get_bins() == {BinAddress(9, 9)}


def test_bin_selection_clear_empties_selection():
    """Verify that clear removes all currently selected bins."""
    sel = BinSelection()
    sel.set_bins([(0, 1), (2, 3)])

    # Verify clear resets the selection to the empty set.
    sel.clear()
    assert sel.get_bins() == set()


def test_bin_selection_get_bins_returns_copy_not_alias():
    """Verify that get_bins returns a defensive copy rather than the internal set."""
    sel = BinSelection()
    sel.set_bins([(1, 1), (2, 2)])

    # Mutate the returned set to verify callers do not receive a live alias.
    bins_copy = sel.get_bins()
    bins_copy.add(BinAddress(999, 999))

    # Verify the internal selection remains unchanged after external mutation.
    assert sel.get_bins() == {BinAddress(1, 1), BinAddress(2, 2)}
    assert BinAddress(999, 999) not in sel.get_bins()


def test_bin_selection_set_bins_accepts_multiple_iterable_input_forms():
    """Verify that set_bins accepts common iterable forms of bin-address input."""
    sel = BinSelection()

    # Verify generator input is consumed and normalized correctly.
    sel.set_bins((i, i + 1) for i in range(3))
    assert sel.get_bins() == {
        BinAddress(0, 1),
        BinAddress(1, 2),
        BinAddress(2, 3),
    }

    # Verify set input is also accepted and normalized correctly.
    sel.set_bins({(5, 6), (7, 8)})
    assert sel.get_bins() == {BinAddress(5, 6), BinAddress(7, 8)}

    # Verify BinAddress instances can be provided directly without conversion issues.
    sel.set_bins([BinAddress(3, 4), BinAddress(3, 4), BinAddress(8, 9)])
    assert sel.get_bins() == {BinAddress(3, 4), BinAddress(8, 9)}
