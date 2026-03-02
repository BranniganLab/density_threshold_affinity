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

from enum import Enum
from dataclasses import dataclass, field


class SelectionOperation(Enum):
    """Enumeration of supported selection operations."""

    REPLACE = "replace"
    ADD = "add"
    SUBTRACT = "subtract"


@dataclass
class SelectorDragState:
    """
    State for tracking a single click/drag gesture.

    Attributes
    ----------
    drag_start
        (r, theta) coordinates of the mouse at the moment the drag began,
        or ``None`` if no drag is in progress.
    last_theta
        Most recently processed (unwrapped) theta for the drag, used to make
        theta continuous across the 0/2π discontinuity.
    current_preview_bins
        The most recently computed preview selection (final selection that would
        be committed if the gesture ended now). ``None`` indicates that no valid
        preview has been computed yet.
    mods
        Modifier set latched at the start of the gesture. Values are strings
        such as ``"shift"`` and ``"control"``. This value should not change
        mid-drag.
    """

    drag_start: tuple[float, float] | None = None
    last_theta: float | None = None
    current_preview_bins: set[tuple[int, int]] | None = None
    mods: frozenset[str] = field(default_factory=frozenset)
