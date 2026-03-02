#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Transient GUI interaction state for Matplotlib-based site selection.

This module defines lightweight, short-lived state objects used by
Matplotlib controllers to track the lifecycle of selection gestures
(mouse press, drag, and release).

The classes here:
- represent user intent (replace/add/remove) derived from modifier keys,
- capture per-gesture drag state and preview geometry,
- exist only during active interactions.

They do not perform rendering or modify persistent selection state.
All drawing is handled by renderers, and all committed selection state
lives in dta.core.selection.

These objects are internal to the Matplotlib GUI layer and are not part
of the public analysis or core APIs.
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
        pass

    def redo(self):
        """Future redo function."""
        pass
