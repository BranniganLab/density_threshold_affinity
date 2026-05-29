#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Per-gesture interaction state for GUI-based selection.

This module defines lightweight state containers used during an active
mouse interaction (press → drag → release) within a `SiteSelector`.
These objects track only the minimal metadata required to interpret a
gesture; they do not perform selection logic and do not store any
preview or committed selection results.

Overview
--------
A selection gesture is initiated by the controller layer
(`SiteSelectorManager`), which determines the active selector and
latches the selection operation (e.g., REPLACE, ADD, SUBTRACT) based
on modifier keys at press time. That latched operation is passed into
the selector and stored for the duration of the gesture.

The classes in this module support that process by maintaining:

- The drag origin in data coordinates
- Angular continuity information across motion events
- The latched selection operation for the active gesture

Contents
--------
SelectionOperation
    Enumeration of selection modes (REPLACE, ADD, SUBTRACT). These
    operations define how a computed bin set should be applied to the
    existing selection.

SelectorDragState
    Container for per-gesture drag state. Tracks the drag origin,
    most recent theta value (for angular continuity), and the latched
    selection operation for the active gesture.

Responsibilities
----------------
- Provide simple, reusable containers for per-gesture state
- Ensure consistent initialization and reset semantics across gestures
- Maintain invariants needed for interpreting motion events

Non-responsibilities
--------------------
- Do not compute or store preview selections
- Do not modify or query the selection model
- Do not interpret keyboard modifiers or determine selection mode

Design Notes
------------
All preview and committed selection state is owned by `SiteSelector`.
This module intentionally remains minimal and policy-free so that
gesture orchestration (controller) and selection semantics (selector)
remain cleanly separated.
"""

from __future__ import annotations

from enum import Enum
from dataclasses import dataclass
from dta.bin_logic.utils import Coordinate


class SelectionOperation(Enum):
    """Enumeration of supported selection operations."""

    REPLACE = "replace"
    ADD = "add"
    SUBTRACT = "subtract"


@dataclass
class SelectorDragState:
    """
    Lightweight container for per-gesture drag state within a SiteSelector.

    This class tracks the minimal transient state required to interpret a
    single mouse drag gesture in polar coordinates. It does not perform any
    selection logic and does not store preview or committed bin sets.

    Responsibilities
    ----------------
    - Store the drag origin (`drag_start`) in data coordinates.
    - Track the most recent theta value (`last_theta`) to maintain angular
      continuity during motion (e.g., across 2π wrap-around).
    - Hold the latched `SelectionOperation` for the duration of the gesture.

    Design Notes
    ------------
    The latched operation is determined by the controller layer
    (`SiteSelectorManager`) at press time and passed into `start_drag(...)`.
    This ensures consistent behavior across the entire gesture, independent
    of subsequent keyboard state changes.

    All preview and committed selection state is owned by `SiteSelector`.
    """

    drag_start: Coordinate | None = None
    last_theta: float | None = None
    _operation: SelectionOperation = SelectionOperation.REPLACE

    @property
    def operation(self) -> SelectionOperation:
        """Get self.operation without exposing setter."""
        return self._operation

    def start_drag(self, at: Coordinate, *, operation: SelectionOperation) -> None:
        """Initialize gesture state at press-time."""
        self.drag_start = Coordinate(*at)
        self.last_theta = self.drag_start.theta_coord
        self._operation = operation

    def reset(self) -> None:
        """Reset to the 'no active gesture' state."""
        self.drag_start = None
        self.last_theta = None
        self._operation = SelectionOperation.REPLACE
