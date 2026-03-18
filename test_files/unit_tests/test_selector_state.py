#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for selection-operation enums and selector drag-state lifecycle."""

import pytest

from dta.bin_logic.utils import Coordinate
from dta.gui.selector_state import SelectionOperation, SelectorDragState


def test_selection_operation_enum_exposes_expected_distinct_values():
    """Verify that SelectionOperation defines the expected distinct string values."""
    # Verify each enum member exposes the public value consumed by callers.
    assert SelectionOperation.REPLACE.value == "replace"
    assert SelectionOperation.ADD.value == "add"
    assert SelectionOperation.SUBTRACT.value == "subtract"

    # Verify the enum members are distinct options rather than aliases.
    assert len(
        {
            SelectionOperation.REPLACE,
            SelectionOperation.ADD,
            SelectionOperation.SUBTRACT,
        }
    ) == 3


def test_selector_drag_state_defaults_are_empty():
    """Verify that a new SelectorDragState starts with no active drag state."""
    s = SelectorDragState()

    # Verify no drag has started yet, so positional tracking is empty.
    assert s.drag_start is None
    assert s.last_theta is None

    # Verify the default semantic selection behavior is replace.
    assert s.operation is SelectionOperation.REPLACE


def test_selector_drag_state_start_drag_latches_operation_and_initializes_position():
    """Verify that start_drag stores the drag origin, seeds theta tracking, and latches the resolved operation."""
    s = SelectorDragState()

    s.start_drag(
        at=(1.25, 0.5),
        operation=SelectionOperation.ADD,
    )

    # Verify the drag origin is coerced to a Coordinate with the expected values.
    assert isinstance(s.drag_start, Coordinate)
    assert s.drag_start == Coordinate(r_coord=1.25, theta_coord=0.5)

    # Verify theta tracking is initialized from the drag start coordinate.
    assert s.last_theta == 0.5

    # Verify the resolved semantic operation is latched exactly as provided.
    assert s.operation is SelectionOperation.ADD


def test_selector_drag_state_reset_clears_active_gesture_state():
    """Verify that reset clears drag position tracking and restores REPLACE mode."""
    s = SelectorDragState()
    s.start_drag(
        at=Coordinate(r_coord=2.0, theta_coord=1.0),
        operation=SelectionOperation.SUBTRACT,
    )

    s.reset()

    assert s.drag_start is None
    assert s.last_theta is None
    assert s.operation is SelectionOperation.REPLACE


def test_selector_drag_state_operation_cannot_be_reassigned_directly():
    """Verify that the latched operation cannot be reassigned directly mid-gesture."""
    s = SelectorDragState()
    s.start_drag(
        at=Coordinate(r_coord=2.0, theta_coord=1.0),
        operation=SelectionOperation.ADD,
    )

    # Verify external callers cannot overwrite the latched operation directly.
    with pytest.raises(AttributeError):
        s.operation = SelectionOperation.SUBTRACT
