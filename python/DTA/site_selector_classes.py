#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 16:22:10 2026.

@author: js2746
"""

from dataclasses import dataclass, field
from typing import Set
import numpy as np
import matplotlib.pyplot as plt
from DTA.utils import unwrap_theta
from DTA.polar_bin_classes import PolarBinGrid, PolarBinRenderer


@dataclass
class SelectorState:
    """
    Save the state for each SiteSelector.

    Parameters
    ----------
    r : tuple[float, float]
        Radial coordinates of the edge endpoints.
    theta : tuple[float, float]
        Angular coordinates of the edge endpoints.
    """

    selected_bins: Set[tuple] = field(default_factory=set)
    selected_artists: list = field(default_factory=list)
    hover_artists: list = field(default_factory=list)
    drag_start: tuple = None
    last_theta: float = None
    mode: str = "replace"


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

    def __init__(self, ax, theta_edges, r_edges, **plotting_args):
        """
        Create a SiteSelector.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Polar axes to attach to.
        theta_edges, r_edges : array-like
            Bin edges.
        plotting_args : dict
            Dictionary of matplotlib-recognized keywords for use with ax.plot.
        """
        self.ax = ax
        self.grid = PolarBinGrid(theta_edges, r_edges)
        self.renderer = PolarBinRenderer(ax)
        self.plotting_args = plotting_args
        self.state = SelectorState()

    def on_activate(self):
        """
        Activate this selector for interaction.

        This method prepares the selector to receive
        mouse events. No selection state is modified.
        """
        self.state.drag_start = None
        self.state.last_theta = None

    def on_deactivate(self):
        """
        Deactivate this selector without destroying state.

        Clears only transient hover and drag state.
        """
        self.state.drag_start = None
        self.state.last_theta = None
        self._clear_artists(self.state.hover_artists)

    def get_selected_bins(self):
        """
        Return selected bins.

        Returns
        -------
        list[tuple[int, int]]
        """
        return sorted(self.state.selected_bins)

    def clear_selection(self):
        """Clear all selected bins and their artists."""
        self.state.selected_bins.clear()
        self._clear_artists(self.state.selected_artists)
        self._clear_artists(self.state.hover_artists)
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
        self._clear_artists(self.state.selected_artists)
        edges = self.grid.exposed_edges(self.state.selected_bins)
        self.state.selected_artists.extend(
            self.renderer.draw_edges(
                edges,
                **self.plotting_args
            )
        )

    def _draw_hover(self, bins):
        """
        Draw hover preview outlines.

        Parameters
        ----------
        bins : iterable of (ri, ti)
        """
        plotting_args = {
            "color": "orange",
            "lw": 1.5,
            "zorder": self.plotting_args["zorder"] + 1,
        }
        self._clear_artists(self.state.hover_artists)
        edges = self.grid.exposed_edges(bins)
        self.state.hover_artists.extend(
            self.renderer.draw_edges(
                edges,
                **plotting_args
            )
        )

    def on_press(self, event):
        """Handle mouse press event."""
        self.state.drag_start = (event.ydata, event.xdata)
        self.state.last_theta = event.xdata

        if event.key == "shift":
            self.state.mode = "add"
        elif event.key == "control":
            self.state.mode = "subtract"
        else:
            self.state.mode = "replace"
            self.state.selected_bins.clear()
            self._clear_artists(self.state.selected_artists)

    def on_motion(self, event):
        """Handle mouse drag event."""
        if self.state.drag_start is None:
            return

        theta = unwrap_theta(self.state.last_theta, event.xdata)
        self.state.last_theta = theta

        bins = self._get_bins_from_drag(
            self.state.drag_start[0], self.state.drag_start[1],
            event.ydata, theta
        )

        bins = set(bins)
        if self.state.mode == "replace":
            temp = bins
        elif self.state.mode == "add":
            temp = self.state.selected_bins | bins
        else:
            temp = self.state.selected_bins - bins

        self._draw_hover(temp)
        self.ax.figure.canvas.draw_idle()

    def on_release(self, event):
        """Handle mouse release event."""
        if self.state.drag_start is None:
            return

        theta = unwrap_theta(self.state.last_theta, event.xdata)

        bins = set(self._get_bins_from_drag(
            self.state.drag_start[0], self.state.drag_start[1],
            event.ydata, theta
        ))

        if self.state.mode in ("replace", "add"):
            self.state.selected_bins.update(bins)
        else:
            self.state.selected_bins.difference_update(bins)

        self._draw_committed()
        self._clear_artists(self.state.hover_artists)

        self.state.drag_start = None
        self.state.last_theta = None
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
            idx = self.grid.map_coord_to_idx(r1, theta1)
            return [idx] if idx else []

        return self.grid.bins_in_region(r0, theta0, r1, theta1)


def example_usage():
    """
    Demonstrate two SiteSelectors on two different Axes objects in same Figure.

    The example creates:
    - A figure with two polar Axes
    - One pcolormesh per Axes
    - One SiteSelector per Axes
    - A shared SiteSelectorManager

    Users can interact with each Axes independently.
    """
    # ------------------------------------------------------------------
    # Create bin edges (20 x 48)
    # ------------------------------------------------------------------
    theta_edges = np.linspace(0, 2 * np.pi, 49)
    r_edges = np.linspace(0.0, 1.0, 21)

    # Meshgrid for pcolormesh (matplotlib wants theta, r)
    theta, r = np.meshgrid(theta_edges, r_edges)

    # Example data for each axes
    z = np.random.rand(20, 48)

    # ------------------------------------------------------------------
    # Create figure and axes
    # ------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(
        1, 2,
        subplot_kw={"projection": "polar"},
        figsize=(10, 5),
        constrained_layout=True,
    )

    # ------------------------------------------------------------------
    # Draw pcolormeshes
    # ------------------------------------------------------------------
    pcm1 = ax1.pcolormesh(
        theta, r, z,
        shading="auto",
        cmap="viridis",
    )
    ax1.set_title("Dataset A")

    pcm2 = ax2.pcolormesh(
        theta, r, z,
        shading="auto",
        cmap="plasma",
    )
    ax2.set_title("Dataset B")

    plotting_args = {
        "color": "red",
        "lw": 2.5,
        "zorder": 20,
    }

    # Optional colorbars
    fig.colorbar(pcm1, ax=ax1, pad=0.1)
    fig.colorbar(pcm2, ax=ax2, pad=0.1)

    # ------------------------------------------------------------------
    # Create SiteSelectors
    # ------------------------------------------------------------------
    selector1 = SiteSelector(
        ax1,
        theta_edges=theta_edges,
        r_edges=r_edges,
        **plotting_args
    )

    selector2 = SiteSelector(
        ax2,
        theta_edges=theta_edges,
        r_edges=r_edges,
        **plotting_args
    )

    # ------------------------------------------------------------------
    # Create and configure the manager
    # ------------------------------------------------------------------
    manager = SiteSelectorManager(fig)

    manager.register(selector1, active=True)
    manager.register(selector2, active=True)

    # ------------------------------------------------------------------
    # Show
    # ------------------------------------------------------------------
    return fig, (ax1, ax2), manager
