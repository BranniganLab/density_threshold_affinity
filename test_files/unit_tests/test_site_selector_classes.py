#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 16:01:53 2026

@author: js2746
"""

from dataclasses import dataclass
import numpy as np

import matplotlib.pyplot as plt

from dta.gui import SiteSelector, SiteSelectorManager


@dataclass
class FakeMouseEvent:
    inaxes: object
    xdata: float | None
    ydata: float | None
    key: str | None = None
    guiEvent: object | None = None
    buttons: tuple[int, ...] | None = (1,)  # emulate "left button down"


def _make_selector(ax):
    theta_edges = np.linspace(0, 2 * np.pi, 9)  # 8 bins
    r_edges = np.linspace(0, 1, 5)            # 4 bins
    return SiteSelector(ax, theta_edges, r_edges, plot_kwargs={"zorder": 20})


def test_selector_freezes_when_mouse_moves_to_other_axes():
    fig, (axA, axB) = plt.subplots(1, 2, subplot_kw={"projection": "polar"})
    selA = _make_selector(axA)

    # Latch replace explicitly
    selA._mods = frozenset()

    selA.on_press(FakeMouseEvent(inaxes=axA, xdata=0.2, ydata=0.2))
    assert selA.drag_tracker.drag_start is not None

    # create a preview
    selA.on_motion(FakeMouseEvent(inaxes=axA, xdata=0.4, ydata=0.4))
    before = selA.drag_tracker.current_preview_bins
    assert before is not None

    # move in Axes B -> should not update preview
    selA.on_motion(FakeMouseEvent(inaxes=axB, xdata=1.2, ydata=0.8))
    assert selA.drag_tracker.current_preview_bins == before


def test_release_commits_last_preview_exactly():
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    # seed committed selection
    sel.selection.set({(0, 0), (1, 1)})

    # latch subtract
    sel._mods = frozenset({"control"})
    sel.on_press(FakeMouseEvent(inaxes=ax, xdata=0.3, ydata=0.3))

    # drag to create preview
    sel.on_motion(FakeMouseEvent(inaxes=ax, xdata=0.6, ydata=0.6))
    preview = sel.drag_tracker.current_preview_bins
    assert preview is not None

    # release outside data coords (xdata/ydata None is fine)
    sel.on_release(FakeMouseEvent(inaxes=None, xdata=None, ydata=None, buttons=None))

    assert sel.selection.snapshot() == frozenset(preview)


def test_add_and_subtract_semantics_across_multiple_drags():
    """
    Coverage item (2): subtract/add preview semantics across multiple drags.
    """
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    sel = _make_selector(ax)

    # Start with some committed bins
    initial = {(0, 0), (0, 1), (1, 1)}
    sel.selection.set(initial)

    # --- SUBTRACT gesture: ctrl drag ---
    sel._mods = frozenset({"control"})
    sel.on_press(FakeMouseEvent(inaxes=ax, xdata=0.2, ydata=0.2))
    sel.on_motion(FakeMouseEvent(inaxes=ax, xdata=0.6, ydata=0.6))
    preview_sub = sel.drag_tracker.current_preview_bins
    assert preview_sub is not None

    sel.on_release(FakeMouseEvent(inaxes=None, xdata=None, ydata=None, buttons=None))
    committed_after_sub = set(sel.selection.get_bins())

    # The committed selection should equal the subtract preview exactly
    assert committed_after_sub == set(preview_sub)

    # --- ADD gesture: shift drag ---
    sel._mods = frozenset({"shift"})
    sel.on_press(FakeMouseEvent(inaxes=ax, xdata=1.0, ydata=0.8))
    sel.on_motion(FakeMouseEvent(inaxes=ax, xdata=1.4, ydata=0.9))
    preview_add = sel.drag_tracker.current_preview_bins
    assert preview_add is not None

    sel.on_release(FakeMouseEvent(inaxes=None, xdata=None, ydata=None, buttons=None))
    committed_after_add = set(sel.selection.get_bins())

    # The committed selection should equal the add preview exactly
    assert committed_after_add == set(preview_add)


def test_integration_freeze_outside_axes_but_commit_on_release_outside():
    """
    Coverage item (3):
    Start drag in A, leave axes (or go into B), do not update preview,
    then release outside and ensure commit matches last hover preview from A.
    """
    fig, (axA, axB) = plt.subplots(1, 2, subplot_kw={"projection": "polar"})
    selA = _make_selector(axA)
    selB = _make_selector(axB)

    mgr = SiteSelectorManager(fig)
    mgr.register(selA, active=True)
    mgr.register(selB, active=True)

    # Start drag on A via manager
    press = FakeMouseEvent(inaxes=axA, xdata=0.2, ydata=0.2)
    mgr._on_press_event(press)

    # Move inside A to generate preview
    mgr._on_motion_event(FakeMouseEvent(inaxes=axA, xdata=0.5, ydata=0.5))
    preview_before = selA.drag_tracker.current_preview_bins
    assert preview_before is not None

    # Move over B — selector A must freeze (no updates)
    mgr._on_motion_event(FakeMouseEvent(inaxes=axB, xdata=1.5, ydata=0.9))
    assert selA.drag_tracker.current_preview_bins == preview_before

    # Release outside axes — manager routes to drag owner
    mgr._on_release_event(FakeMouseEvent(inaxes=None, xdata=None, ydata=None, buttons=None))
    assert selA.drag_tracker.drag_start is None

    # Commit should equal last preview from A
    assert selA.selection.snapshot() == frozenset(preview_before)
