#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 16:22:10 2026.

@author: js2746
"""
import numpy as np
import matplotlib.pyplot as plt
from enum import Enum
from DTA.utils import unwrap_theta
from DTA.polar_bin_classes import PolarBinGrid, PolarBinRenderer, BinSelectionModel


class SelectionOperation(Enum):
    """Enumeration of supported selection operations."""

    REPLACE = "replace"
    ADD = "add"
    SUBTRACT = "subtract"


class SiteSelector:
    """
    Interactive controller for selecting polar bins on a single Axes.

    This class converts mouse gestures into semantic selection operations,
    updates a selection model, and renders both hover previews and
    committed selections.
    """

    def __init__(self, ax, theta_edges, r_edges, *, color="red", lw=2, zorder=10):
        """
        Create a SiteSelector object.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Polar axes used for interaction and drawing.
        theta_edges, r_edges : array-like
            Bin edge definitions.
        color : str
            Color used to draw committed selections.
        lw : float
            Line width used to draw committed selections.
        zorder : int
            Drawing order for committed selections.
        """
        self.ax = ax
        self.grid = PolarBinGrid(theta_edges, r_edges)
        self.renderer = PolarBinRenderer(ax)
        self.model = BinSelectionModel()

        self.color = color
        self.lw = lw
        self.zorder = zorder

        self._hover_artists = []
        self._selected_artists = []

        self._drag_start = None
        self._last_theta = None
        self._operation = SelectionOperation.REPLACE

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def on_activate(self):
        """
        Prepare the selector to receive user input.

        This does not modify the current selection.
        """
        self._drag_start = None
        self._last_theta = None

    def on_deactivate(self):
        """Disable interaction without modifying the selection."""
        self._clear_artists(self._hover_artists)
        self._drag_start = None
        self._last_theta = None

    # ------------------------------------------------------------------
    # Event handlers
    # ------------------------------------------------------------------

    def on_press(self, event):
        """
        Begin a selection gesture.

        Modifier keys determine how the selection will be applied:
        - no modifier: replace
        - Shift: add
        - Control: subtract
        """
        self._drag_start = (event.ydata, event.xdata)
        self._last_theta = event.xdata

        if event.key == "shift":
            self._operation = SelectionOperation.ADD
        elif event.key == "control":
            self._operation = SelectionOperation.SUBTRACT
        else:
            self._operation = SelectionOperation.REPLACE
            self._clear_artists(self._selected_artists)
            self.model.clear()

    def on_motion(self, event):
        """Update the hover preview while dragging."""
        if self._drag_start is None:
            return

        theta = unwrap_theta(self._last_theta, event.xdata)
        self._last_theta = theta

        bins = self._bins_from_drag(
            self._drag_start, (event.ydata, theta)
        )

        preview_bins = self._apply_preview(bins)
        self._draw_hover(preview_bins)

        self.ax.figure.canvas.draw_idle()

    def on_release(self, event):
        """Finalize a selection gesture and commit the result."""
        if self._drag_start is None:
            return

        before = self.model.snapshot()

        theta = unwrap_theta(self._last_theta, event.xdata)
        bins = self._bins_from_drag(
            self._drag_start, (event.ydata, theta)
        )

        self._apply_commit(bins)
        after = self.model.snapshot()

        self.on_selection_committed(before, after)

        self._draw_committed()
        self._clear_artists(self._hover_artists)

        self._drag_start = None
        self._last_theta = None
        self.ax.figure.canvas.draw_idle()

    # ------------------------------------------------------------------
    # Selection logic
    # ------------------------------------------------------------------

    def _bins_from_drag(self, start, end):
        """
        Determine which bins are covered by a click or drag gesture.

        Parameters
        ----------
        start : tuple[float, float]
            (r, theta) of the drag start.
        end : tuple[float, float]
            (r, theta) of the drag end.

        Returns
        -------
        set[tuple[int, int]]
            Selected bin indices.
        """
        r0, t0 = start
        r1, t1 = end

        if abs(r1 - r0) < 1e-8 and abs(t1 - t0) < 1e-8:
            idx = self.grid.bin_at(r1, t1)
            return {idx} if idx is not None else set()

        return set(self.grid.bins_in_region(r0, t0, r1, t1))

    def _apply_preview(self, bins):
        """Compute the preview selection without mutating the model."""
        current = self.model.bins()

        if self._operation is SelectionOperation.REPLACE:
            return bins
        if self._operation is SelectionOperation.ADD:
            return current | bins
        if self._operation is SelectionOperation.SUBTRACT:
            return current - bins

        return current

    def _apply_commit(self, bins):
        """Apply the selection operation to the model."""
        if self._operation is SelectionOperation.REPLACE:
            self.model.set(bins)
        elif self._operation is SelectionOperation.ADD:
            self.model.add(bins)
        elif self._operation is SelectionOperation.SUBTRACT:
            self.model.remove(bins)

    # ------------------------------------------------------------------
    # Rendering
    # ------------------------------------------------------------------

    def _draw_hover(self, bins):
        """Draw a temporary hover preview for a set of bins."""
        self._clear_artists(self._hover_artists)
        edges = self.grid.exposed_edges(bins)

        self._hover_artists.extend(
            self.renderer.draw_edges(
                edges,
                color="orange",
                lw=1.5,
                zorder=self.zorder + 1,
            )
        )

    def _draw_committed(self):
        """Draw the committed selection."""
        self._clear_artists(self._selected_artists)
        edges = self.grid.exposed_edges(self.model.bins())

        self._selected_artists.extend(
            self.renderer.draw_edges(
                edges,
                color=self.color,
                lw=self.lw,
                zorder=self.zorder,
            )
        )

    def _clear_artists(self, artists):
        """Remove matplotlib artists from the Axes and clear the list."""
        for artist in artists:
            artist.remove()
        artists.clear()

    # ------------------------------------------------------------------
    # Undo hook
    # ------------------------------------------------------------------

    def on_selection_committed(self, before, after):
        """
        Save space for future implementation of 'undo' logic.

        Parameters
        ----------
        before : frozenset
            Selection state before the operation.
        after : frozenset
            Selection state after the operation.
        """
        pass


class SiteSelectorManager:
    """
    Event router that manages multiple SiteSelectors within a Figure.

    Exactly one selector per Axes is active at a time.
    """

    def __init__(self, fig):
        """
        Create a SiteSelectorManager and tie it to a Matplotlib Figure.

        Parameters
        ----------
        fig : matplotlib.figure.Figure
        """
        self.fig = fig
        self._selectors = {}
        self._active = {}

        self._cids = [
            fig.canvas.mpl_connect("button_press_event", self._dispatch("on_press")),
            fig.canvas.mpl_connect("motion_notify_event", self._dispatch("on_motion")),
            fig.canvas.mpl_connect("button_release_event", self._dispatch("on_release")),
        ]

    def register(self, selector, *, active=False):
        """
        Register a selector with the manager.

        Parameters
        ----------
        selector : SiteSelector
        active : bool
            Whether this selector should initially receive input.
        """
        ax = selector.ax
        self._selectors.setdefault(ax, []).append(selector)
        if active or ax not in self._active:
            self.set_active(selector)

    def set_active(self, selector):
        """Activate a selector for its Axes."""
        ax = selector.ax
        current = self._active.get(ax)
        if current is selector:
            return

        if current:
            current.on_deactivate()

        self._active[ax] = selector
        selector.on_activate()

    def _dispatch(self, method):
        """Create an event handler that forwards events to the active selector."""
        def handler(event):
            selector = self._active.get(event.inaxes)
            if selector:
                getattr(selector, method)(event)
        return handler


def example_usage():
    """
    Create a figure with two polar pcolormesh plots and attach one SiteSelector to each Axes.

    Each SiteSelector operates independently and allows interactive
    selection of polar bins on its associated Axes.
    """
    # ------------------------------------------------------------------
    # Define polar bin edges
    # ------------------------------------------------------------------
    n_theta = 48
    n_r = 20

    theta_edges = np.linspace(0.0, 2.0 * np.pi, n_theta + 1)
    r_edges = np.linspace(0.0, 1.0, n_r + 1)

    # Meshgrid required by pcolormesh (theta, r ordering)
    Theta, R = np.meshgrid(theta_edges, r_edges)

    # Example data for each plot
    data_a = np.random.rand(n_r, n_theta)
    data_b = np.random.rand(n_r, n_theta)

    # ------------------------------------------------------------------
    # Create figure and Axes
    # ------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(
        1, 2,
        subplot_kw={"projection": "polar"},
        figsize=(10, 5),
        constrained_layout=True,
    )

    # ------------------------------------------------------------------
    # Draw pcolormesh plots
    # ------------------------------------------------------------------
    pcm1 = ax1.pcolormesh(
        Theta,
        R,
        data_a,
        shading="auto",
        cmap="viridis",
    )
    ax1.set_title("Polar Data A")
    fig.colorbar(pcm1, ax=ax1, pad=0.1)

    pcm2 = ax2.pcolormesh(
        Theta,
        R,
        data_b,
        shading="auto",
        cmap="plasma",
    )
    ax2.set_title("Polar Data B")
    fig.colorbar(pcm2, ax=ax2, pad=0.1)

    # ------------------------------------------------------------------
    # Create SiteSelectors (one per Axes)
    # ------------------------------------------------------------------
    selector_a = SiteSelector(
        ax1,
        theta_edges=theta_edges,
        r_edges=r_edges,
        color="red",
        lw=2.0,
        zorder=20,
    )

    selector_b = SiteSelector(
        ax2,
        theta_edges=theta_edges,
        r_edges=r_edges,
        color="cyan",
        lw=2.0,
        zorder=20,
    )

    # ------------------------------------------------------------------
    # Register selectors with a manager
    # ------------------------------------------------------------------
    manager = SiteSelectorManager(fig)
    manager.register(selector_a, active=True)
    manager.register(selector_b, active=True)

    return fig, (ax1, ax2), manager
