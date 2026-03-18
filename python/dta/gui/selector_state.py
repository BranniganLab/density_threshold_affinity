#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ephemeral interaction state for Matplotlib-driven selection gestures.

This module contains small, data-centric types used by the Matplotlib GUI layer
to track a *single* selection gesture (press → drag → release). These objects
are intentionally short-lived and UI-facing: they capture latched modifier keys,
drag continuity bookkeeping, and the current preview result during an active
interaction.

Contents
--------
SelectionOperation
    Semantic operation requested by the user (replace/add/subtract), typically
    derived from modifier keys at press-time.

SelectorDragState
    Per-gesture state latched and updated during interaction:
    - drag start location in (r, theta)
    - theta unwrapping continuity (last_theta)
    - the most recent preview bin set that would be committed on release
    - modifier set frozen at gesture start

Non-goals
---------
- No rendering: these types do not touch Matplotlib artists/canvases.
- No persistence: committed selection state is owned by the selection model in
  the core/bin-logic layer (e.g., dta.bin_logic), not here.
- No domain validation: address validity and invariants belong to the domain
  model and grid/geometry code.

This module is an internal implementation detail of the Matplotlib GUI and is
not part of the public analysis API.
"""

from __future__ import annotations

from enum import Enum
from dataclasses import dataclass
from dta.bin_logic import Coordinate


class SelectionOperation(Enum):
    """Enumeration of supported selection operations."""

    REPLACE = "replace"
    ADD = "add"
    SUBTRACT = "subtract"


@dataclass
class SelectorDragState:
    """
    State for tracking a single click/drag gesture.

    drag_start is stored as a Coordinate once a drag begins.
    operation is the semantic selection action resolved by the controller
    at gesture start.
    """

    drag_start: Coordinate | None = None
    last_theta: float | None = None
    operation: SelectionOperation = SelectionOperation.REPLACE

    def start_drag(self, at: Coordinate, *, operation: SelectionOperation) -> None:
        """Initialize gesture state at press-time."""
        self.drag_start = Coordinate(*at)
        self.last_theta = self.drag_start.theta_coord
        self.operation = operation

    def update_theta(self, theta_unwrapped: float) -> None:
        """Update continuity bookkeeping."""
        self.last_theta = theta_unwrapped

    def clear(self) -> None:
        """Reset to the 'no active gesture' state."""
        self.drag_start = None
        self.last_theta = None
        self.operation = SelectionOperation.REPLACE
