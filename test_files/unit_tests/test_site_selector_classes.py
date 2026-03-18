#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for SiteSelector and SiteSelectorManager."""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np

from dta.gui import SiteSelector, SiteSelectorManager
from dta.gui.selector_state import SelectionOperation


@dataclass
class FakeMouseEvent:
    """Minimal mouse-event stub for selector/controller unit tests."""

    inaxes: object
    xdata: float | None
    ydata: float | None
    key: str | None = None
    guiEvent: object | None = None
    buttons: tuple[int, ...] | None = (1,)


def _make_selector(ax):
    """Construct a selector over a small polar grid for test use."""
    theta_edges = np.linspace(0, 2 * np.pi, 9)  # 8 bins
    r_edges = np.linspace(0, 1, 5)              # 4 bins
    return SiteSelector(ax, theta_edges, r_edges, plot_kwargs={"zorder": 20})


def test_selector_freezes_when_mouse_moves_to_other_axes():
    """Preview should stop updating when an active drag leaves the selector axes."""
    fig, (ax_a, ax_b) = plt.subplots(1, 2, subplot_kw={"projection": "polar"})
    sel_a = _make_selector(ax_a)

    sel_a.on_press(
        FakeMouseEvent(inaxes=ax_a, xdata=0.2, ydata=0.2),
        SelectionOperation.REPLACE,
    )
    assert sel_a.drag_tracker.drag_start is not None
    assert sel_a.drag_tracker.operation is SelectionOperation.REPLACE

    # Establish a preview while the cursor remains inside the owning axes.
    sel_a.on_motion(FakeMouseEvent(inaxes=ax_a, xdata=0.4, ydata=0.4))
    before = sel_a.current_preview_bins
    assert before is not None

    # Moving into a different axes should freeze the preview.
    sel_a.on_motion(FakeMouseEvent(inaxes=ax_b, xdata=1.2, ydata=0.8))
    assert sel_a.current_preview_bins == before


def test_release_commits_last_preview_exactly():
    """Release should commit exactly the most recent preview selection."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    # Seed an existing committed selection.
    sel.selection.set_bins({(0, 0), (1, 1)})

    # Start a subtract gesture and build a preview.
    sel.on_press(
        FakeMouseEvent(inaxes=ax, xdata=0.3, ydata=0.3),
        SelectionOperation.SUBTRACT,
    )
    assert sel.drag_tracker.operation is SelectionOperation.SUBTRACT

    sel.on_motion(FakeMouseEvent(inaxes=ax, xdata=0.6, ydata=0.6))
    preview = sel.current_preview_bins
    assert preview is not None

    # Release outside the axes; release does not need xdata/ydata.
    sel.on_release(FakeMouseEvent(inaxes=None, xdata=None, ydata=None, buttons=None))

    assert sel.selection.get_bins() == preview
    assert sel.drag_tracker.drag_start is None
    assert sel.drag_tracker.last_theta is None
    assert sel.drag_tracker.operation is SelectionOperation.REPLACE
    assert sel.current_preview_bins is None


def test_add_and_subtract_semantics_across_multiple_drags():
    """Subtract and add gestures should each commit exactly the preview they compute."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    # Start with some committed bins already selected.
    initial = {(0, 0), (0, 1), (1, 1)}
    sel.selection.set_bins(initial)

    # First gesture: subtract.
    sel.on_press(
        FakeMouseEvent(inaxes=ax, xdata=0.2, ydata=0.2),
        SelectionOperation.SUBTRACT,
    )
    assert sel.drag_tracker.operation is SelectionOperation.SUBTRACT

    sel.on_motion(FakeMouseEvent(inaxes=ax, xdata=0.6, ydata=0.6))
    preview_sub = sel.current_preview_bins
    assert preview_sub is not None

    sel.on_release(FakeMouseEvent(inaxes=None, xdata=None, ydata=None, buttons=None))
    committed_after_sub = set(sel.selection.get_bins())

    assert committed_after_sub == set(preview_sub)
    assert sel.current_preview_bins is None
    assert sel.drag_tracker.operation is SelectionOperation.REPLACE

    # Second gesture: add.
    sel.on_press(
        FakeMouseEvent(inaxes=ax, xdata=1.0, ydata=0.8),
        SelectionOperation.ADD,
    )
    assert sel.drag_tracker.operation is SelectionOperation.ADD

    sel.on_motion(FakeMouseEvent(inaxes=ax, xdata=1.4, ydata=0.9))
    preview_add = sel.current_preview_bins
    assert preview_add is not None

    sel.on_release(FakeMouseEvent(inaxes=None, xdata=None, ydata=None, buttons=None))
    committed_after_add = set(sel.selection.get_bins())

    assert committed_after_add == set(preview_add)
    assert sel.current_preview_bins is None
    assert sel.drag_tracker.operation is SelectionOperation.REPLACE


def test_manager_latches_replace_when_no_modifiers():
    """A plain press should latch REPLACE mode on the drag owner."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    mgr = SiteSelectorManager(fig)
    mgr.register(sel, active=True)

    mgr._on_press_event(FakeMouseEvent(inaxes=ax, xdata=0.2, ydata=0.2))

    assert mgr._drag_owner is sel
    assert sel.drag_tracker.operation is SelectionOperation.REPLACE


def test_manager_latches_add_when_shift_is_pressed():
    """Shift should latch ADD mode for the duration of the gesture."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    mgr = SiteSelectorManager(fig)
    mgr.register(sel, active=True)

    mgr._on_press_event(
        FakeMouseEvent(inaxes=ax, xdata=0.2, ydata=0.2, key="shift")
    )

    assert mgr._drag_owner is sel
    assert sel.drag_tracker.operation is SelectionOperation.ADD


def test_manager_latches_subtract_when_control_is_pressed():
    """Control should latch SUBTRACT mode for the duration of the gesture."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    mgr = SiteSelectorManager(fig)
    mgr.register(sel, active=True)

    mgr._on_press_event(
        FakeMouseEvent(inaxes=ax, xdata=0.2, ydata=0.2, key="control")
    )

    assert mgr._drag_owner is sel
    assert sel.drag_tracker.operation is SelectionOperation.SUBTRACT


def test_manager_prefers_shift_when_shift_and_control_are_both_present():
    """If both modifiers are present, shift semantics should win."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    mgr = SiteSelectorManager(fig)
    mgr.register(sel, active=True)

    mgr._on_press_event(
        FakeMouseEvent(inaxes=ax, xdata=0.2, ydata=0.2, key="shift+control")
    )

    assert mgr._drag_owner is sel
    assert sel.drag_tracker.operation is SelectionOperation.ADD


def test_manager_reads_gui_event_modifier_flags():
    """Manager should honor modifier flags supplied through guiEvent."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    mgr = SiteSelectorManager(fig)
    mgr.register(sel, active=True)

    mgr._on_press_event(
        FakeMouseEvent(
            inaxes=ax,
            xdata=0.2,
            ydata=0.2,
            guiEvent={"ctrlKey": True},
        )
    )

    assert mgr._drag_owner is sel
    assert sel.drag_tracker.operation is SelectionOperation.SUBTRACT


def test_integration_freeze_outside_axes_but_commit_on_release_outside():
    """Manager should keep routing to the drag owner until release."""
    fig, (ax_a, ax_b) = plt.subplots(1, 2, subplot_kw={"projection": "polar"})
    sel_a = _make_selector(ax_a)
    sel_b = _make_selector(ax_b)

    mgr = SiteSelectorManager(fig)
    mgr.register(sel_a, active=True)
    mgr.register(sel_b, active=True)

    # Start drag on A via the manager.
    mgr._on_press_event(FakeMouseEvent(inaxes=ax_a, xdata=0.2, ydata=0.2))
    assert mgr._drag_owner is sel_a
    assert sel_a.drag_tracker.operation is SelectionOperation.REPLACE

    # Move inside A to create the preview that should later be committed.
    mgr._on_motion_event(FakeMouseEvent(inaxes=ax_a, xdata=0.5, ydata=0.5))
    preview_before = sel_a.current_preview_bins
    assert preview_before is not None

    # Move over B; A should freeze and retain its last preview.
    mgr._on_motion_event(FakeMouseEvent(inaxes=ax_b, xdata=1.5, ydata=0.9))
    assert sel_a.current_preview_bins == preview_before

    # Release outside all axes; the manager should still route to A.
    mgr._on_release_event(
        FakeMouseEvent(inaxes=None, xdata=None, ydata=None, buttons=None)
    )

    assert mgr._drag_owner is None
    assert sel_a.drag_tracker.drag_start is None
    assert sel_a.drag_tracker.last_theta is None
    assert sel_a.drag_tracker.operation is SelectionOperation.REPLACE
    assert sel_a.current_preview_bins is None
    assert sel_a.selection.get_bins() == preview_before
