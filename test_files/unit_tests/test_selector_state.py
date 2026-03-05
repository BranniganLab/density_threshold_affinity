#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 11:27:51 2026

@author: js2746
"""

import pytest

from dta.bin_logic import Coordinate
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
    assert s.mods == frozenset()
    assert isinstance(s.mods, frozenset)


def test_selector_drag_state_start_drag_latches_mods_and_coerces_coordinate():
    s = SelectorDragState()

    s.start_drag(
        at=(1.25, 0.5),
        mods=frozenset({"shift", "control"}),
    )

    assert isinstance(s.drag_start, Coordinate)
    assert s.drag_start == Coordinate(r_coord=1.25, theta_coord=0.5)

    # start_drag initializes last_theta from drag_start.theta
    assert s.last_theta == 0.5
    assert s.mods == frozenset({"shift", "control"})

    # mods is a frozenset, so it should be immutable
    with pytest.raises(AttributeError):
        s.mods.add("alt")  # type: ignore[attr-defined]


def test_selector_drag_state_update_theta_and_clear():
    s = SelectorDragState()
    s.start_drag(at=Coordinate(r_coord=2.0, theta_coord=1.0), mods=frozenset({"shift"}))

    s.update_theta(6.1)
    assert s.last_theta == 6.1

    s.clear()
    assert s.drag_start is None
    assert s.last_theta is None
    assert s.mods == frozenset()
