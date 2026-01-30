#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 16:22:10 2026.

@author: js2746
"""
from enum import Enum
from dataclasses import dataclass, field
import numpy as np
import matplotlib.pyplot as plt
from DTA.utils import unwrap_theta
from DTA.polar_bin_classes import PolarBinGrid, PolarBinRenderer, BinSelectionModel


class SelectionOperation(Enum):
    """Enumeration of supported selection operations."""

    REPLACE = "replace"
    ADD = "add"
    SUBTRACT = "subtract"


@dataclass
class SelectorDragState:
    """State for tracking drag events."""

    drag_start: tuple[float, float] | None = None
    last_theta: float | None = None
    last_preview_bins: set[tuple[int, int]] | None = None


@dataclass
class SelectorDrawState:
    """State for tracking what is drawn."""

    selected_artists: list = field(default_factory=list)
    hover_artists: list = field(default_factory=list)


class SiteSelector:
    """
    Interactive controller for selecting polar bins on a single Axes.

    This class converts mouse gestures into semantic selection operations,
    updates a selection model, and renders both hover previews and
    committed selections.
    """

    def __init__(self, ax, theta_edges, r_edges, plot_kwargs={}):
        """
        Create a SiteSelector object.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Polar axes used for interaction and drawing.
        theta_edges, r_edges : array-like
            Bin edge definitions.
        plot_kwargs : dict, OPTIONAL
            Dictionary of matplotlib plotting keywords for drawing bin edges.
        """
        self.ax = ax
        self.grid = PolarBinGrid(theta_edges, r_edges)
        self.renderer = PolarBinRenderer(ax, plot_kwargs)
        self.model = BinSelectionModel()
        self.draw_tracker = SelectorDrawState()
        self.drag_tracker = SelectorDragState()
        self.operation = SelectionOperation.REPLACE

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def on_activate(self):
        """
        Prepare the selector to receive user input.

        This does not modify the current selection.
        """
        self.drag_tracker.drag_start = None
        self.drag_tracker.last_theta = None
        self.drag_tracker.last_preview_bins = None

    def on_deactivate(self):
        """Disable interaction without modifying the selection."""
        self._clear_artists(self.draw_tracker.hover_artists)
        self.drag_tracker.drag_start = None
        self.drag_tracker.last_theta = None
        self.drag_tracker.last_preview_bins = None

    # ------------------------------------------------------------------
    # Event handlers
    # ------------------------------------------------------------------

    def on_press(self, event):
        """Begin a selection gesture."""
        if event.inaxes is not self.ax:
            return

        self.drag_tracker.drag_start = (event.ydata, event.xdata)
        self.drag_tracker.last_theta = event.xdata
        self.drag_tracker.last_preview_bins = None

        mods = getattr(self, "_mods", frozenset())

        if "shift" in mods:
            self.operation = SelectionOperation.ADD
        elif "control" in mods:
            self.operation = SelectionOperation.SUBTRACT
        else:
            self.operation = SelectionOperation.REPLACE

        # Optional: prime a click-preview immediately if xdata/ydata are valid
        if event.xdata is not None and event.ydata is not None:
            start = (event.ydata, event.xdata)
            bins = self._bins_from_drag(start, start)
            preview_bins = self._apply_preview(bins)
            self.drag_tracker.last_preview_bins = preview_bins
            self._draw_hover(preview_bins)
            self.ax.figure.canvas.draw_idle()

    def on_motion(self, event):
        """Update the hover preview while dragging."""
        if self.drag_tracker.drag_start is None:
            return

        # If we're outside the axes (or outside data), do not update the hover.
        if event.inaxes is not self.ax:
            return
        if event.xdata is None or event.ydata is None:
            return

        theta = unwrap_theta(self.drag_tracker.last_theta, event.xdata)
        self.drag_tracker.last_theta = theta

        bins = self._bins_from_drag(
            self.drag_tracker.drag_start, (event.ydata, theta)
        )

        preview_bins = self._apply_preview(bins)

        # Cache exactly what we drew, so release can commit the same set.
        self.drag_tracker.last_preview_bins = preview_bins

        self._draw_hover(preview_bins)
        self.ax.figure.canvas.draw_idle()

    def on_release(self, event):
        """Finalize a selection gesture and commit the result."""
        if self.drag_tracker.drag_start is None:
            return

        before = self.model.snapshot()

        preview_bins = self.drag_tracker.last_preview_bins

        if preview_bins is None:
            # No hover was ever drawn (e.g., click with no motion or first motion invalid).
            # Fall back to click-at-press using the press coordinates.
            r0, t0 = self.drag_tracker.drag_start
            idx = self.grid.bin_at(r0, t0)
            bins = {idx} if idx is not None else set()

            # Here bins is a delta for the current operation, so use the delta-commit path.
            self._apply_commit(bins)
        else:
            # preview_bins is already the FINAL desired selection. Commit it as-is.
            self._commit_preview_selection(preview_bins)

        after = self.model.snapshot()
        self.on_selection_committed(before, after)

        self._draw_committed()
        self._clear_artists(self.draw_tracker.hover_artists)

        self.drag_tracker.drag_start = None
        self.drag_tracker.last_theta = None
        self.drag_tracker.last_preview_bins = None
        self._mods = frozenset()
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
            # treat as single click, not drag event
            idx = self.grid.bin_at(r1, t1)
            return {idx} if idx is not None else set()

        # else: treat as drag event
        return set(self.grid.bins_in_region(r0, t0, r1, t1))

    def _apply_preview(self, bins):
        """Compute the preview selection without mutating the model."""
        current = self.model.bins()

        if self.operation is SelectionOperation.REPLACE:
            return bins
        if self.operation is SelectionOperation.ADD:
            return current | bins
        if self.operation is SelectionOperation.SUBTRACT:
            return current - bins

        return current

    def _apply_commit(self, bins):
        """Apply the selection operation to the model."""
        if self.operation is SelectionOperation.REPLACE:
            self.model.set(bins)
        elif self.operation is SelectionOperation.ADD:
            self.model.add(bins)
        elif self.operation is SelectionOperation.SUBTRACT:
            self.model.remove(bins)

    def _commit_preview_selection(self, preview_bins: set[tuple[int, int]]):
        """
        Commit a *final selection set* (not a delta).

        This is required when on_release commits the last hover preview, because
        preview_bins already includes the effect of ADD/SUBTRACT/REPLACE.
        """
        self.model.set(preview_bins)

    # ------------------------------------------------------------------
    # Rendering
    # ------------------------------------------------------------------

    def _draw_hover(self, bins):
        """Draw a temporary hover preview for a set of bins."""
        self._clear_artists(self.draw_tracker.hover_artists)
        edges = self.grid.exposed_edges(bins)

        hover_kwargs = {
            "color": "orange",
            "lw": 1.5,
            "zorder": self.renderer.plot_kwargs['zorder'] + 1
        }

        self.draw_tracker.hover_artists.extend(
            self.renderer.draw_edges(
                edges,
                hover_kwargs
            )
        )

    def _draw_committed(self):
        """Draw the committed selection."""
        self._clear_artists(self.draw_tracker.selected_artists)
        edges = self.grid.exposed_edges(self.model.bins())

        self.draw_tracker.selected_artists.extend(
            self.renderer.draw_edges(
                edges
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


class SiteSelectorManager:
    """
    Event router that manages multiple SiteSelectors within a Figure.

    Exactly one selector per Axes is active at a time.

    ipympl / widget backend note:
    - MouseEvent.key is unreliable.
    - key_press/key_release are unreliable.
    - The most reliable modifier info is in event.guiEvent (often dict-like).
    """

    def __init__(self, fig):
        self.fig = fig
        self._selectors = {}
        self._active = {}
        self._drag_owner = None

        self._cids = [
            fig.canvas.mpl_connect("button_press_event", self._dispatch("on_press")),
            fig.canvas.mpl_connect("motion_notify_event", self._dispatch("on_motion")),
            fig.canvas.mpl_connect("button_release_event", self._dispatch("on_release")),
        ]

    def register(self, selector, *, active=False):
        ax = selector.ax
        self._selectors.setdefault(ax, []).append(selector)
        if active or ax not in self._active:
            self.set_active(selector)

    def set_active(self, selector):
        ax = selector.ax
        current = self._active.get(ax)
        if current is selector:
            return
        if current:
            current.on_deactivate()
        self._active[ax] = selector
        selector.on_activate()

    def _mods_from_mouse_event(self, event) -> set[str]:
        """
        Extract modifiers from the mouse event itself.

        In ipympl, event.guiEvent may be a dict but sometimes lacks modifier fields.
        Only return early if we actually found something in guiEvent; otherwise
        fall back to event.key.
        """
        mods: set[str] = set()
        ge = getattr(event, "guiEvent", None)

        def get_bool(obj, key: str) -> bool:
            if isinstance(obj, dict):
                return bool(obj.get(key, False))
            return bool(getattr(obj, key, False))

        # Try guiEvent first (works only if it actually includes modifier flags)
        if ge is not None:
            try:
                if get_bool(ge, "shiftKey"):
                    mods.add("shift")
                if get_bool(ge, "ctrlKey") or get_bool(ge, "metaKey"):
                    mods.add("control")
                # IMPORTANT: do NOT return unless we found something
                if mods:
                    return mods
            except Exception:
                # If guiEvent isn't shaped as expected, fall back below
                pass

        # Fallback: event.key string (often the only thing available in widget backend)
        k = (getattr(event, "key", None) or "").lower()
        if "shift" in k:
            mods.add("shift")
        if "control" in k or "ctrl" in k or "meta" in k:
            mods.add("control")

        return mods

    def _dispatch(self, method):
        def handler(event):
            # If a drag is in progress, always route motion/release to drag owner.
            if self._drag_owner is not None and method in ("on_motion", "on_release"):
                getattr(self._drag_owner, method)(event)
                if method == "on_release":
                    # reset latch after the gesture ends
                    self._drag_owner._mods = frozenset()
                    self._drag_owner = None
                return

            selector = self._active.get(getattr(event, "inaxes", None))
            if not selector:
                return

            if method == "on_press":
                self._drag_owner = selector
                mods = self._mods_from_mouse_event(event)
                selector._mods = frozenset(mods)

            getattr(selector, method)(event)

        return handler


def example_usage():
    """
    Create a figure with two polar pcolormesh plots and attach one SiteSelector to each Axes.

    Each SiteSelector operates independently and allows interactive
    selection of polar bins on its associated Axes.
    """
    # ------------------------------------------------------------------
    # Define polar bin edges (20x48 lattice bins)
    # ------------------------------------------------------------------
    theta_edges = np.linspace(0.0, 2.0 * np.pi, 49)
    r_edges = np.linspace(0.0, 1.0, 21)

    # Meshgrid required by pcolormesh (theta, r ordering)
    theta, r = np.meshgrid(theta_edges, r_edges)

    # Example data for each plot
    data = np.random.rand(20, 48)

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
        theta,
        r,
        data,
        shading="auto",
        cmap="viridis",
    )
    ax1.set_title("Polar Data A")
    fig.colorbar(pcm1, ax=ax1, pad=0.1)

    pcm2 = ax2.pcolormesh(
        theta,
        r,
        data,
        shading="auto",
        cmap="plasma",
    )
    ax2.set_title("Polar Data B")
    fig.colorbar(pcm2, ax=ax2, pad=0.1)

    # ------------------------------------------------------------------
    # Create SiteSelectors (one per Axes)
    # ------------------------------------------------------------------
    plotting_kwargs = {
        "color": "red",
        "lw": 2.0,
        "zorder": 20
    }

    selector_a = SiteSelector(
        ax1,
        theta_edges=theta_edges,
        r_edges=r_edges,
        plot_kwargs=plotting_kwargs,
    )

    plotting_kwargs['color'] = 'cyan'

    selector_b = SiteSelector(
        ax2,
        theta_edges=theta_edges,
        r_edges=r_edges,
        plot_kwargs=plotting_kwargs,
    )

    # ------------------------------------------------------------------
    # Register selectors with a manager
    # ------------------------------------------------------------------
    manager = SiteSelectorManager(fig)
    manager.register(selector_a, active=True)
    manager.register(selector_b, active=True)

    return fig, (ax1, ax2), manager
