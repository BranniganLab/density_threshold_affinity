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


@dataclass
class ObjectGuiEvent:
    """Object-style guiEvent stub exposing modifier attributes."""

    shiftKey: bool = False
    ctrlKey: bool = False
    metaKey: bool = False


class DummyArtist:
    """Simple stand-in for a Matplotlib artist that records removal."""

    def __init__(self):
        self.removed = False

    def remove(self):
        """Mark this artist as removed."""
        self.removed = True


def _make_selector(ax):
    """Construct a selector over a small polar grid for test use."""
    theta_edges = np.linspace(0, 2 * np.pi, 9)  # 8 angular bins
    r_edges = np.linspace(0, 1, 5)              # 4 radial bins
    return SiteSelector(ax, theta_edges, r_edges, plot_kwargs={"zorder": 20})


# ============================================================================
# SiteSelector direct tests
# ============================================================================


def test_site_selector_on_press_returns_false_outside_axes():
    """Presses outside the selector axes should be ignored."""
    fig, (ax_a, ax_b) = plt.subplots(1, 2, subplot_kw={"projection": "polar"})
    sel = _make_selector(ax_a)

    updated = sel.on_press(
        FakeMouseEvent(inaxes=ax_b, xdata=0.2, ydata=0.2),
        SelectionOperation.REPLACE,
    )

    assert updated is False
    assert sel.drag_tracker.drag_start is None
    assert sel.current_preview_bins is None


def test_site_selector_on_press_returns_false_for_missing_data_coords():
    """Presses with missing data coordinates should be ignored."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    updated = sel.on_press(
        FakeMouseEvent(inaxes=ax, xdata=None, ydata=0.2),
        SelectionOperation.REPLACE,
    )

    assert updated is False
    assert sel.drag_tracker.drag_start is None
    assert sel.current_preview_bins is None


def test_site_selector_on_press_latches_operation_and_creates_initial_preview():
    """Valid press should initialize drag state and establish an initial preview."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    updated = sel.on_press(
        FakeMouseEvent(inaxes=ax, xdata=0.2, ydata=0.3),
        SelectionOperation.ADD,
    )

    assert updated is True
    assert sel.drag_tracker.drag_start is not None
    assert sel.drag_tracker.operation is SelectionOperation.ADD
    assert sel.current_preview_bins is not None


def test_site_selector_on_motion_returns_false_without_active_drag():
    """Motion events should be ignored when no drag gesture is active."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    updated = sel.on_motion(FakeMouseEvent(inaxes=ax, xdata=0.4, ydata=0.4))

    assert updated is False
    assert sel.current_preview_bins is None


def test_site_selector_on_motion_returns_false_outside_axes_and_freezes_preview():
    """Motion outside the owning axes should not update the existing preview."""
    fig, (ax_a, ax_b) = plt.subplots(1, 2, subplot_kw={"projection": "polar"})
    sel = _make_selector(ax_a)

    sel.on_press(
        FakeMouseEvent(inaxes=ax_a, xdata=0.2, ydata=0.2),
        SelectionOperation.REPLACE,
    )
    sel.on_motion(FakeMouseEvent(inaxes=ax_a, xdata=0.4, ydata=0.4))
    preview_before = sel.current_preview_bins

    updated = sel.on_motion(FakeMouseEvent(inaxes=ax_b, xdata=0.6, ydata=0.6))

    assert updated is False
    assert sel.current_preview_bins == preview_before


def test_site_selector_on_motion_returns_false_for_missing_data_coords():
    """Motion with missing data coordinates should be ignored during a drag."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    sel.on_press(
        FakeMouseEvent(inaxes=ax, xdata=0.2, ydata=0.2),
        SelectionOperation.REPLACE,
    )
    preview_before = sel.current_preview_bins

    updated = sel.on_motion(FakeMouseEvent(inaxes=ax, xdata=None, ydata=0.4))

    assert updated is False
    assert sel.current_preview_bins == preview_before


def test_site_selector_on_motion_updates_preview_and_last_theta():
    """Valid motion should update both the preview bins and angular tracking."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    sel.on_press(
        FakeMouseEvent(inaxes=ax, xdata=0.2, ydata=0.2),
        SelectionOperation.REPLACE,
    )
    last_theta_before = sel.drag_tracker.last_theta

    updated = sel.on_motion(FakeMouseEvent(inaxes=ax, xdata=0.6, ydata=0.6))

    assert updated is True
    assert sel.current_preview_bins is not None
    assert sel.drag_tracker.last_theta != last_theta_before


def test_site_selector_on_release_returns_false_without_active_drag():
    """Release should be ignored when no drag gesture is active."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    updated = sel.on_release(
        FakeMouseEvent(inaxes=None, xdata=None, ydata=None, buttons=None)
    )

    assert updated is False
    assert sel.drag_tracker.drag_start is None
    assert sel.current_preview_bins is None


def test_site_selector_release_commits_last_preview_exactly():
    """Release should commit exactly the most recent preview selection."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    sel.selection.set_bins({(0, 0), (1, 1)})

    sel.on_press(
        FakeMouseEvent(inaxes=ax, xdata=0.3, ydata=0.3),
        SelectionOperation.SUBTRACT,
    )
    sel.on_motion(FakeMouseEvent(inaxes=ax, xdata=0.6, ydata=0.6))
    preview = sel.current_preview_bins

    updated = sel.on_release(
        FakeMouseEvent(inaxes=None, xdata=None, ydata=None, buttons=None)
    )

    assert updated is True
    assert preview is not None
    assert sel.selection.get_bins() == preview
    assert sel.drag_tracker.drag_start is None
    assert sel.drag_tracker.last_theta is None
    assert sel.drag_tracker.operation is SelectionOperation.REPLACE
    assert sel.current_preview_bins is None


def test_site_selector_add_and_subtract_semantics_across_multiple_drags():
    """Subtract and add gestures should each commit exactly the preview they compute."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    initial = {(0, 0), (0, 1), (1, 1)}
    sel.selection.set_bins(initial)

    sel.on_press(
        FakeMouseEvent(inaxes=ax, xdata=0.2, ydata=0.2),
        SelectionOperation.SUBTRACT,
    )
    sel.on_motion(FakeMouseEvent(inaxes=ax, xdata=0.6, ydata=0.6))
    preview_sub = sel.current_preview_bins

    sel.on_release(FakeMouseEvent(inaxes=None, xdata=None, ydata=None, buttons=None))
    committed_after_sub = set(sel.selection.get_bins())

    assert preview_sub is not None
    assert committed_after_sub == set(preview_sub)
    assert sel.current_preview_bins is None
    assert sel.drag_tracker.operation is SelectionOperation.REPLACE

    sel.on_press(
        FakeMouseEvent(inaxes=ax, xdata=1.0, ydata=0.8),
        SelectionOperation.ADD,
    )
    sel.on_motion(FakeMouseEvent(inaxes=ax, xdata=1.4, ydata=0.9))
    preview_add = sel.current_preview_bins

    sel.on_release(FakeMouseEvent(inaxes=None, xdata=None, ydata=None, buttons=None))
    committed_after_add = set(sel.selection.get_bins())

    assert preview_add is not None
    assert committed_after_add == set(preview_add)
    assert sel.current_preview_bins is None
    assert sel.drag_tracker.operation is SelectionOperation.REPLACE


def test_site_selector_calculate_preview_bins_replace_mode():
    """REPLACE preview should equal the drag-implied bin set."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    sel.selection.set_bins({(0, 0), (1, 1)})
    sel.drag_tracker.start_drag(
        sel.grid.bin_center(0, 0),
        operation=SelectionOperation.REPLACE,
    )

    bins = {(2, 2), (3, 1)}
    preview = sel._calculate_preview_bins(bins)

    assert preview == bins


def test_site_selector_calculate_preview_bins_add_mode():
    """ADD preview should union the drag-implied bins with the current selection."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    current = {(0, 0), (1, 1)}
    bins = {(1, 1), (2, 2)}
    sel.selection.set_bins(current)
    sel.drag_tracker.start_drag(
        sel.grid.bin_center(0, 0),
        operation=SelectionOperation.ADD,
    )

    preview = sel._calculate_preview_bins(bins)

    assert preview == current | bins


def test_site_selector_calculate_preview_bins_subtract_mode():
    """SUBTRACT preview should remove drag-implied bins from the current selection."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    current = {(0, 0), (1, 1), (2, 2)}
    bins = {(1, 1), (3, 3)}
    sel.selection.set_bins(current)
    sel.drag_tracker.start_drag(
        sel.grid.bin_center(0, 0),
        operation=SelectionOperation.SUBTRACT,
    )

    preview = sel._calculate_preview_bins(bins)

    assert preview == {(0, 0), (2, 2)}


def test_site_selector_on_activate_resets_transient_state():
    """Activation should clear drag metadata and preview state only."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    sel.selection.set_bins({(0, 0)})
    sel.on_press(
        FakeMouseEvent(inaxes=ax, xdata=0.2, ydata=0.2),
        SelectionOperation.ADD,
    )
    assert sel.drag_tracker.drag_start is not None
    assert sel.current_preview_bins is not None

    sel.on_activate()

    assert sel.drag_tracker.drag_start is None
    assert sel.drag_tracker.last_theta is None
    assert sel.drag_tracker.operation is SelectionOperation.REPLACE
    assert sel.current_preview_bins is None
    assert sel.selection.get_bins() == {(0, 0)}


def test_site_selector_on_deactivate_clears_hover_artists_and_transient_state():
    """Deactivation should remove hover artists and reset transient gesture state."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    artist_a = DummyArtist()
    artist_b = DummyArtist()
    sel.renderer.hover_artists.extend([artist_a, artist_b])

    sel.on_press(
        FakeMouseEvent(inaxes=ax, xdata=0.2, ydata=0.2),
        SelectionOperation.REPLACE,
    )
    assert sel.drag_tracker.drag_start is not None

    sel.on_deactivate()

    assert artist_a.removed is True
    assert artist_b.removed is True
    assert sel.renderer.hover_artists == []
    assert sel.drag_tracker.drag_start is None
    assert sel.drag_tracker.last_theta is None
    assert sel.drag_tracker.operation is SelectionOperation.REPLACE
    assert sel.current_preview_bins is None


# ============================================================================
# SiteSelectorManager direct tests
# ============================================================================


def test_site_selector_manager_set_active_deactivates_previous_selector():
    """Switching the active selector on one axes should deactivate the previous one."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel_a = _make_selector(ax)
    sel_b = _make_selector(ax)

    manager = SiteSelectorManager(fig)

    sel_a.on_press(
        FakeMouseEvent(inaxes=ax, xdata=0.2, ydata=0.2),
        SelectionOperation.ADD,
    )
    sel_a.on_motion(FakeMouseEvent(inaxes=ax, xdata=0.4, ydata=0.4))
    assert sel_a.drag_tracker.drag_start is not None
    assert sel_a.current_preview_bins is not None

    manager.register(sel_a, active=True)
    manager.set_active(sel_b)

    assert manager._active[ax] is sel_b
    assert sel_a.drag_tracker.drag_start is None
    assert sel_a.current_preview_bins is None


def test_site_selector_manager_mods_from_dict_gui_event_shift():
    """Modifier extraction should recognize shift from dict-style guiEvent."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    manager = SiteSelectorManager(fig)

    mods = manager._mods_from_mouse_event(
        FakeMouseEvent(
            inaxes=ax,
            xdata=0.2,
            ydata=0.2,
            guiEvent={"shiftKey": True},
        )
    )

    assert mods == {"shift"}


def test_site_selector_manager_mods_from_dict_gui_event_ctrl_and_meta():
    """Modifier extraction should treat ctrl or meta as control semantics."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    manager = SiteSelectorManager(fig)

    mods_ctrl = manager._mods_from_mouse_event(
        FakeMouseEvent(
            inaxes=ax,
            xdata=0.2,
            ydata=0.2,
            guiEvent={"ctrlKey": True},
        )
    )
    mods_meta = manager._mods_from_mouse_event(
        FakeMouseEvent(
            inaxes=ax,
            xdata=0.2,
            ydata=0.2,
            guiEvent={"metaKey": True},
        )
    )

    assert mods_ctrl == {"control"}
    assert mods_meta == {"control"}


def test_site_selector_manager_mods_from_object_gui_event():
    """Modifier extraction should recognize object-style guiEvent attributes."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    manager = SiteSelectorManager(fig)

    mods = manager._mods_from_mouse_event(
        FakeMouseEvent(
            inaxes=ax,
            xdata=0.2,
            ydata=0.2,
            guiEvent=ObjectGuiEvent(shiftKey=True, ctrlKey=True),
        )
    )

    assert mods == {"shift", "control"}


def test_site_selector_manager_mods_from_key_fallback():
    """Modifier extraction should fall back to parsing event.key strings."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    manager = SiteSelectorManager(fig)

    mods = manager._mods_from_mouse_event(
        FakeMouseEvent(
            inaxes=ax,
            xdata=0.2,
            ydata=0.2,
            key="shift+ctrl",
        )
    )

    assert mods == {"shift", "control"}


def test_site_selector_manager_latches_replace_when_no_modifiers():
    """A plain press should latch REPLACE mode on the drag owner."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    manager = SiteSelectorManager(fig)
    manager.register(sel, active=True)

    manager._on_press_event(FakeMouseEvent(inaxes=ax, xdata=0.2, ydata=0.2))

    assert manager._drag_owner is sel
    assert sel.drag_tracker.operation is SelectionOperation.REPLACE


def test_site_selector_manager_latches_add_when_shift_is_pressed():
    """Shift should latch ADD mode for the duration of the gesture."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    manager = SiteSelectorManager(fig)
    manager.register(sel, active=True)

    manager._on_press_event(
        FakeMouseEvent(inaxes=ax, xdata=0.2, ydata=0.2, key="shift")
    )

    assert manager._drag_owner is sel
    assert sel.drag_tracker.operation is SelectionOperation.ADD


def test_site_selector_manager_latches_subtract_when_control_is_pressed():
    """Control should latch SUBTRACT mode for the duration of the gesture."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    manager = SiteSelectorManager(fig)
    manager.register(sel, active=True)

    manager._on_press_event(
        FakeMouseEvent(inaxes=ax, xdata=0.2, ydata=0.2, key="control")
    )

    assert manager._drag_owner is sel
    assert sel.drag_tracker.operation is SelectionOperation.SUBTRACT


def test_site_selector_manager_prefers_shift_when_shift_and_control_are_both_present():
    """If both modifiers are present, shift semantics should win."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    manager = SiteSelectorManager(fig)
    manager.register(sel, active=True)

    manager._on_press_event(
        FakeMouseEvent(inaxes=ax, xdata=0.2, ydata=0.2, key="shift+control")
    )

    assert manager._drag_owner is sel
    assert sel.drag_tracker.operation is SelectionOperation.ADD


def test_site_selector_manager_motion_returns_early_without_drag_owner():
    """Motion callback should do nothing when no selector owns the drag."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    manager = SiteSelectorManager(fig)

    manager._on_motion_event(FakeMouseEvent(inaxes=ax, xdata=0.2, ydata=0.2))

    assert manager._drag_owner is None


def test_site_selector_manager_release_returns_early_without_drag_owner():
    """Release callback should do nothing when no selector owns the drag."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    manager = SiteSelectorManager(fig)

    manager._on_release_event(
        FakeMouseEvent(inaxes=ax, xdata=0.2, ydata=0.2, buttons=None)
    )

    assert manager._drag_owner is None


# ============================================================================
# End-to-end interaction tests
# ============================================================================


def test_integration_freeze_outside_axes_but_commit_on_release_outside():
    """Manager should keep routing to the drag owner until release."""
    fig, (ax_a, ax_b) = plt.subplots(1, 2, subplot_kw={"projection": "polar"})
    sel_a = _make_selector(ax_a)
    sel_b = _make_selector(ax_b)

    manager = SiteSelectorManager(fig)
    manager.register(sel_a, active=True)
    manager.register(sel_b, active=True)

    manager._on_press_event(FakeMouseEvent(inaxes=ax_a, xdata=0.2, ydata=0.2))
    assert manager._drag_owner is sel_a
    assert sel_a.drag_tracker.operation is SelectionOperation.REPLACE

    manager._on_motion_event(FakeMouseEvent(inaxes=ax_a, xdata=0.5, ydata=0.5))
    preview_before = sel_a.current_preview_bins
    assert preview_before is not None

    manager._on_motion_event(FakeMouseEvent(inaxes=ax_b, xdata=1.5, ydata=0.9))
    assert sel_a.current_preview_bins == preview_before

    manager._on_release_event(
        FakeMouseEvent(inaxes=None, xdata=None, ydata=None, buttons=None)
    )

    assert manager._drag_owner is None
    assert sel_a.drag_tracker.drag_start is None
    assert sel_a.drag_tracker.last_theta is None
    assert sel_a.drag_tracker.operation is SelectionOperation.REPLACE
    assert sel_a.current_preview_bins is None
    assert sel_a.selection.get_bins() == preview_before
