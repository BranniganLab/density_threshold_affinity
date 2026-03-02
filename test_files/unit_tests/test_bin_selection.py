#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 14:13:46 2026

@author: js2746
"""

from dta.bin_logic import BinSelection


def test_bin_selection_initially_empty():
    sel = BinSelection()
    assert sel.get_bins() == set()


def test_bin_selection_update_replaces_selection():
    sel = BinSelection()

    sel.set_bins([(0, 1), (2, 3), (0, 1)])  # duplicates should be collapsed
    assert sel.get_bins() == {(0, 1), (2, 3)}

    sel.set_bins([(9, 9)])
    assert sel.get_bins() == {(9, 9)}


def test_bin_selection_clear_empties_selection():
    sel = BinSelection()
    sel.set_bins([(0, 1), (2, 3)])
    sel.clear()
    assert sel.get_bins() == set()


def test_bin_selection_get_bins_returns_copy_not_alias():
    sel = BinSelection()
    sel.set_bins([(1, 1), (2, 2)])

    bins_copy = sel.get_bins()
    bins_copy.add((999, 999))  # mutate the returned set

    # internal state should be unchanged
    assert sel.get_bins() == {(1, 1), (2, 2)}
    assert (999, 999) not in sel.get_bins()


def test_bin_selection_update_accepts_any_iterable():
    sel = BinSelection()

    # generator input
    sel.set_bins((i, i + 1) for i in range(3))
    assert sel.get_bins() == {(0, 1), (1, 2), (2, 3)}

    # set input
    sel.set_bins({(5, 6), (7, 8)})
    assert sel.get_bins() == {(5, 6), (7, 8)}
