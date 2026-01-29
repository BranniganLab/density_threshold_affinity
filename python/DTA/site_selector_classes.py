#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 16:22:10 2026.

@author: js2746
"""
from DTA.utils import unwrap_theta
from DTA.polar_bin_classes import PolarBinGrid, PolarBinRenderer


class SiteSelectorManager:
    """
    Central event router for SiteSelectors.

    One manager should exist per Figure.
    """

    def __init__(self, fig):
        """
        Create a SiteSelectorManager.

        Parameters
        ----------
        fig : matplotlib.figure.Figure
        """
        self.fig = fig
        self._selectors_by_axes = {}
        self._active = {}

        self._connections = [
            fig.canvas.mpl_connect("button_press_event", self._on_press),
            fig.canvas.mpl_connect("motion_notify_event", self._on_motion),
            fig.canvas.mpl_connect("button_release_event", self._on_release),
        ]

    def register(self, selector, *, active=False):
        """
        Register a selector with the manager.

        Parameters
        ----------
        selector : SiteSelector
        active : bool
            Whether this selector should be active initially.
        """
        ax = selector.ax
        self._selectors_by_axes.setdefault(ax, []).append(selector)

        if active or ax not in self._active:
            self.set_active(selector)

    def set_active(self, selector):
        """Set the active selector for its Axes."""
        ax = selector.ax
        current = self._active.get(ax)

        if current is selector:
            return

        if current is not None:
            current.on_deactivate()

        self._active[ax] = selector
        selector.on_activate()

    def _dispatch(self, method, event):
        """Dispatch an event to the active selector."""
        if event.inaxes is None:
            return

        selector = self._active.get(event.inaxes)
        if selector is not None:
            getattr(selector, method)(event)

    def _on_press(self, event):
        """Handle button press events."""
        self._dispatch("on_press", event)

    def _on_motion(self, event):
        """Handle mouse motion events."""
        self._dispatch("on_motion", event)

    def _on_release(self, event):
        """Handle button release events."""
        self._dispatch("on_release", event)

    def disconnect(self):
        """Disconnect all event handlers and clear state."""
        for cid in self._connections:
            self.fig.canvas.mpl_disconnect(cid)
        self._selectors_by_axes.clear()
        self._active.clear()


class SiteSelector:
    """
    Interactive selection layer for polar bins.

    Each SiteSelector:
    - Owns its own selection state
    - Draws persistent outlines on an Axes
    - Can be activated or deactivated without losing state
    """

    def __init__(self, ax, theta_edges, r_edges, *, color="red", lw=2, zorder=10):
        """
        Create a SiteSelector.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Polar axes to attach to.
        theta_edges, r_edges : array-like
            Bin edges.
        color : str
            Color for committed selections.
        lw : float
            Line width for committed selections.
        zorder : int
            Z-order for committed selections.
        """
        self.ax = ax
        self.grid = PolarBinGrid(theta_edges, r_edges)
        self.renderer = PolarBinRenderer(ax)

        self.color = color
        self.lw = lw
        self.zorder = zorder

        self.selected_bins = set()
        self.selected_artists = []
        self.hover_artists = []

        self.drag_start = None
        self._last_theta = None
        self.mode = "replace"

    def on_activate(self):
        """
        Activate this selector for interaction.

        This method prepares the selector to receive
        mouse events. No selection state is modified.
        """
        self.drag_start = None
        self._last_theta = None

    def on_deactivate(self):
        """
        Deactivate this selector without destroying state.

        Clears only transient hover and drag state.
        """
        self.drag_start = None
        self._last_theta = None
        self._clear_artists(self.hover_artists)

    def get_selected_bins(self):
        """
        Return selected bins.

        Returns
        -------
        list[tuple[int, int]]
        """
        return sorted(self.selected_bins)

    def clear_selection(self):
        """Clear all selected bins and their artists."""
        self.selected_bins.clear()
        self._clear_artists(self.selected_artists)
        self._clear_artists(self.hover_artists)
        self.ax.figure.canvas.draw_idle()

    def _clear_artists(self, artists):
        """
        Remove matplotlib artists from the Axes.

        Parameters
        ----------
        artists : list
        """
        for a in artists:
            a.remove()
        artists.clear()

    def _draw_committed(self):
        """Draw committed (confirmed) selection outlines."""
        self._clear_artists(self.selected_artists)
        edges = self.grid.exposed_edges(self.selected_bins)
        self.selected_artists.extend(
            self.renderer.draw_edges(
                edges,
                color=self.color,
                lw=self.lw,
                zorder=self.zorder,
            )
        )

    def _draw_hover(self, bins):
        """
        Draw hover preview outlines.

        Parameters
        ----------
        bins : iterable of (ri, ti)
        """
        self._clear_artists(self.hover_artists)
        edges = self.grid.exposed_edges(bins)
        self.hover_artists.extend(
            self.renderer.draw_edges(
                edges,
                color="orange",
                lw=1.5,
                zorder=self.zorder + 1,
            )
        )

    def on_press(self, event):
        """Handle mouse press event."""
        self.drag_start = (event.ydata, event.xdata)
        self._last_theta = event.xdata

        if event.key == "shift":
            self.mode = "add"
        elif event.key == "control":
            self.mode = "subtract"
        else:
            self.mode = "replace"
            self.selected_bins.clear()
            self._clear_artists(self.selected_artists)

    def on_motion(self, event):
        """Handle mouse drag event."""
        if self.drag_start is None:
            return

        theta = unwrap_theta(self._last_theta, event.xdata)
        self._last_theta = theta

        bins = self._get_bins_from_drag(
            self.drag_start[0], self.drag_start[1],
            event.ydata, theta
        )

        bins = set(bins)
        if self.mode == "replace":
            temp = bins
        elif self.mode == "add":
            temp = self.selected_bins | bins
        else:
            temp = self.selected_bins - bins

        self._draw_hover(temp)
        self.ax.figure.canvas.draw_idle()

    def on_release(self, event):
        """Handle mouse release event."""
        if self.drag_start is None:
            return

        theta = unwrap_theta(self._last_theta, event.xdata)

        bins = set(self._get_bins_from_drag(
            self.drag_start[0], self.drag_start[1],
            event.ydata, theta
        ))

        if self.mode in ("replace", "add"):
            self.selected_bins.update(bins)
        else:
            self.selected_bins.difference_update(bins)

        self._draw_committed()
        self._clear_artists(self.hover_artists)

        self.drag_start = None
        self._last_theta = None
        self.ax.figure.canvas.draw_idle()

    def _get_bins_from_drag(self, r0, theta0, r1, theta1):
        """
        Determine bins from a click or drag.

        Parameters
        ----------
        r0, theta0 : float
            Drag start coordinates.
        r1, theta1 : float
            Drag end coordinates.

        Returns
        -------
        list[tuple[int, int]]
        """
        if abs(theta1 - theta0) < 1e-8 and abs(r1 - r0) < 1e-8:
            # treat this as a single click, not a drag
            idx = self.grid.bin_at(r1, theta1)
            return [idx] if idx else []

        return self.grid.bins_in_region(r0, theta0, r1, theta1)
