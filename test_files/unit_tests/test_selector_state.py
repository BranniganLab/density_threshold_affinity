#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 11:27:51 2026

@author: js2746
"""

import pytest

from dta.gui.selector_state import BinSelection, SelectionOperation, SelectorDragState


def test_selection_operation_enum_values():
    assert SelectionOperation.REPLACE.value == "replace"
    assert SelectionOperation.ADD.value == "add"
    assert SelectionOperation.SUBTRACT.value == "subtract"

    # sanity: enum members are distinct
    assert len({SelectionOperation.REPLACE, SelectionOperation.ADD, SelectionOperation.SUBTRACT}) == 3


def test_selector_drag_state_defaults_are_empty_and_none():
    s = SelectorDragState()

    assert s.drag_start is None
    assert s.last_theta is None
    assert s.current_preview_bins is None
    assert s.mods == frozenset()
    assert isinstance(s.mods, frozenset)


def test_selector_drag_state_allows_setting_fields_and_latches_mods():
    s = SelectorDragState(
        drag_start=(1.25, 0.5),
        last_theta=6.1,
        current_preview_bins={(0, 0), (1, 2)},
        mods=frozenset({"shift", "control"}),
    )

    assert s.drag_start == (1.25, 0.5)
    assert s.last_theta == 6.1
    assert s.current_preview_bins == {(0, 0), (1, 2)}
    assert s.mods == frozenset({"shift", "control"})

    # mods is a frozenset, so it should be immutable
    with pytest.raises(AttributeError):
        s.mods.add("alt")  # type: ignore[attr-defined]


def test_bin_selection_initially_empty():
    sel = BinSelection()
    assert sel.get_bins() == set()


def test_bin_selection_update_replaces_selection():
    sel = BinSelection()

    sel.update([(0, 1), (2, 3), (0, 1)])  # duplicates should be collapsed
    assert sel.get_bins() == {(0, 1), (2, 3)}

    sel.update([(9, 9)])
    assert sel.get_bins() == {(9, 9)}


def test_bin_selection_clear_empties_selection():
    sel = BinSelection()
    sel.update([(0, 1), (2, 3)])
    sel.clear()
    assert sel.get_bins() == set()


def test_bin_selection_get_bins_returns_copy_not_alias():
    sel = BinSelection()
    sel.update([(1, 1), (2, 2)])

    bins_copy = sel.get_bins()
    bins_copy.add((999, 999))  # mutate the returned set

    # internal state should be unchanged
    assert sel.get_bins() == {(1, 1), (2, 2)}
    assert (999, 999) not in sel.get_bins()


def test_bin_selection_update_accepts_any_iterable():
    sel = BinSelection()

    # generator input
    sel.update((i, i + 1) for i in range(3))
    assert sel.get_bins() == {(0, 1), (1, 2), (2, 3)}

    # set input
    sel.update({(5, 6), (7, 8)})
    assert sel.get_bins() == {(5, 6), (7, 8)}
