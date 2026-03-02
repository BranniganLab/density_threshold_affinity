#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 11:27:51 2026

@author: js2746
"""

import pytest

from dta.gui.selector_state import SelectionOperation, SelectorDragState


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
