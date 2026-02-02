#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
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
    last_preview_bins
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
    last_preview_bins: set[tuple[int, int]] | None = None
    mods: frozenset[str] = field(default_factory=frozenset)


@dataclass
class SelectorDrawState:
    """
    State tracking Matplotlib artists created by the selector.

    Attributes
    ----------
    selected_artists
        Artists representing the committed selection.
    hover_artists
        Artists representing the transient hover preview.
    """

    selected_artists: list = field(default_factory=list)
    hover_artists: list = field(default_factory=list)
