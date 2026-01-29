#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 16 13:02:14 2026.

@author: js2746
"""
import numpy as np
import matplotlib.pyplot as plt
from DTA.utils import bin_in_theta_arc, unwrap_theta


# pylint: disable=too-many-instance-attributes
class SiteSelector:
    """
    Interactive bin selector for polar plots.

    Usage:
        fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
        ax.pcolormesh(theta_edges, r_edges, data)
        selector = PolarBinSelector(ax, theta_edges, r_edges)
        plt.show()

        # Access selected bins
        selected = selector.get_selected_bins()
    """

    def __init__(self, ax, theta_edges, r_edges):
        """
        Initialize the polar bin selector.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Polar axes to attach selector to
        theta_edges : array-like
            Angular bin edges (radians)
        r_edges : array-like
            Radial bin edges
        """
        self.ax = ax
        self.theta_edges = np.asarray(theta_edges)
        self.r_edges = np.asarray(r_edges)

        # State
        self.selected_bins = set()
        self.hover_artists = []
        self.selected_artists = []
        self.drag_start = None
        self._last_theta = None
        self.mode = "replace"

        # Connect events
        self.fig = ax.figure
        self._connections = [
            self.fig.canvas.mpl_connect("button_press_event", self._on_press),
            self.fig.canvas.mpl_connect("motion_notify_event", self._on_motion),
            self.fig.canvas.mpl_connect("button_release_event", self._on_release)
        ]

        # Remove matplotlib default content
        self.fig.canvas.toolbar_visible = False
        self.fig.canvas.header_visible = False
        self.fig.canvas.footer_visible = False
        self.fig.canvas.resizable = False

    def get_selected_bins(self):
        """
        Get the currently selected bins.

        Returns
        -------
        list of tuple
            List of (theta_idx, r_idx) tuples
        """
        return sorted(list(self.selected_bins))

    def clear_selection(self):
        """Clear all selected bins and remove artists."""
        self.selected_bins.clear()
        self._clear_artists(self.selected_artists)
        self._clear_artists(self.hover_artists)
        self.fig.canvas.draw_idle()

    def _clear_artists(self, artists):
        """Remove artists from plot and clear list."""
        for a in artists:
            a.remove()
        artists.clear()

    def _bins_in_region(self, theta_start, theta_end, r_start, r_end):
        r_start, r_end = sorted([r_start, r_end])

        bins = []
        for ti in range(len(self.theta_edges) - 1):
            t_low = self.theta_edges[ti]
            t_high = self.theta_edges[ti + 1]
            if bin_in_theta_arc(theta_start, theta_end, t_low, t_high):
                for ri in range(len(self.r_edges) - 1):
                    r_low = self.r_edges[ri]
                    r_high = self.r_edges[ri + 1]
                    if r_high >= r_start and r_low <= r_end:
                        bins.append((ti, ri))
        return bins

    def _draw_outer_edges(self, mask_bins, **plotting_args):
        """
        Draw an outline around the outer boundary of a set of selected polar bins.

        Only edges that border an unselected neighbor (or the plot boundary)
        are drawn.
        """
        if not mask_bins:
            return []

        mask = self._build_bin_mask(mask_bins)
        artists = []

        for ri, ti in zip(*np.where(mask)):
            artists.extend(
                self._draw_bin_edges(mask, ri, ti, **plotting_args)
            )

        return artists

    def _build_bin_mask(self, mask_bins):
        """Return a boolean (r, theta) mask of selected bins."""
        n_r = len(self.r_edges) - 1
        n_t = len(self.theta_edges) - 1
        mask = np.zeros((n_r, n_t), dtype=bool)

        for ti, ri in mask_bins:
            mask[ri, ti] = True

        return mask

    def _draw_bin_edges(self, mask, ri, ti, **plotting_args):
        """
        Draw the visible (outer) edges of a single selected bin.

        Each bin has four possible edges:

                top
             ┌─────────┐
             │         │
        left │ (ri,ti) │ right
             │         │
             └─────────┘
               bottom

        An edge is drawn if:
            - The neighboring bin in that direction is unselected
            - OR the bin lies on the outer boundary of the grid

        This ensures:
            - Continuous outlines
            - No duplicated internal edges

        Parameters
        ----------
        mask : ndarray
            2D boolean array of selected bins.
        ri : int
            r index of bin.
        ti : int
            Theta index of bin.
        color : str
            Line color to be drawn.
        lw : float
            Line width to be drawn.

        Returns
        -------
        list of matplotlib artists to draw.
        """
        edges = self._bin_edges(ri, ti)
        draw = self._determine_if_edges_exposed(mask, ri, ti)

        artists = []

        if draw["top"]:
            artists.append(self._plot_edge(edges["top"], **plotting_args))
        if draw["bottom"]:
            artists.append(self._plot_edge(edges["bottom"], **plotting_args))
        if draw["left"]:
            artists.append(self._plot_edge(edges["left"], **plotting_args))
        if draw["right"]:
            artists.append(self._plot_edge(edges["right"], **plotting_args))

        return artists

    def _bin_edges(self, ri, ti):
        """
        Return the polar line segments corresponding to a bin's edges.

        Parameters
        ----------
        ri : int
            r index of bin.
        ti : int
            Theta index of bin.

        Returns
        -------
        Dictionary of line segments, one for each edge of the bin.
        """
        t0 = self.theta_edges[ti]
        t1 = self.theta_edges[(ti + 1) % (len(self.theta_edges) - 1)]
        r0 = self.r_edges[ri]
        r1 = self.r_edges[ri + 1]

        return {
            "top":    ([t0, t1], [r1, r1]),
            "bottom": ([t0, t1], [r0, r0]),
            "left":   ([t0, t0], [r0, r1]),
            "right":  ([t1, t1], [r0, r1]),
        }

    def _determine_if_edges_exposed(self, mask, ri, ti):
        """
        Determine which edges of a bin are exposed (i.e., neighbor not selected).

        Parameters
        ----------
        mask : ndarray
            2D boolean array of selected bins.
        ri : int
            r index of bin.
        ti : int
            Theta index of bin.

        Returns
        -------
        Dictionary of booleans, one for each edge of the bin.
        """
        n_r, n_t = mask.shape

        return {
            "top":    ri == n_r - 1 or not mask[ri + 1, ti],
            "bottom": ri == 0 or not mask[ri - 1, ti],
            "left":   not mask[ri, (ti - 1) % n_t],
            "right":  not mask[ri, (ti + 1) % n_t],
        }

    def _plot_edge(self, line, **plotting_args):
        """Plot a single polar edge line."""
        theta_vals, r_vals = line
        return self.ax.plot(theta_vals, r_vals, **plotting_args, zorder=10)[0]

    def _get_bins_from_drag(self, theta_start, r_start, theta_end, r_end):
        """
        Determine whether click or drag and return selected bins.

        Parameters
        ----------
        theta_start : float
            Theta coordinate of drag/click start.
        r_start : float
            r coordinate of drag/click start.
        theta_end : float
            Theta coordinate of drag/click end.
        r_end : float
            r coordinate of drag/click end.

        Returns
        -------
        list of tuple
            List of (theta_index, r_index) tuples

        """
        dt = theta_end - theta_start
        dr = r_end - r_start
        threshold = 1e-8

        if abs(dt) < threshold and abs(dr) < threshold:
            # Single bin click, not drag event
            ti = np.searchsorted(self.theta_edges, theta_end % (2 * np.pi), side="right") - 1
            ri = np.searchsorted(self.r_edges, r_end, side="right") - 1
            if 0 <= ti < len(self.theta_edges) - 1 and 0 <= ri < len(self.r_edges) - 1:
                # In-bounds click, return current bin
                return [(ti, ri)]
            # else: out-of-bounds click, return empty list
            return []
        # else: this was a drag event, return all the bins in the drag region
        return self._bins_in_region(theta_start, theta_end, r_start, r_end)

    def _update_hover_display(self, bins):
        """Update drawn region as mouse drags."""
        plot_args = {
            'color': 'orange',
            'lw': 1.5,
        }
        self._clear_artists(self.hover_artists)

        if self.mode == "replace":
            temp_mask = set(bins)
        elif self.mode == "add":
            temp_mask = self.selected_bins.union(bins)
        elif self.mode == "subtract":
            temp_mask = self.selected_bins.difference(bins)
        else:
            temp_mask = self.selected_bins.union(bins)

        self.hover_artists.extend(
            self._draw_outer_edges(temp_mask, **plot_args)
        )
        self.fig.canvas.draw_idle()

    def _update_selected_display(self):
        """Update confirmed selection display."""
        plot_args = {
            'color': 'red',
            'lw': 2,
        }
        self._clear_artists(self.selected_artists)
        self.selected_artists.extend(
            self._draw_outer_edges(self.selected_bins, **plot_args)
        )
        self.fig.canvas.draw_idle()

    def _on_press(self, event):
        """Handle mouse press event."""
        if event.inaxes != self.ax or event.button != 1:
            return

        self.drag_start = (event.xdata, event.ydata)
        self._last_theta = event.xdata

        if event.key == "shift":
            self.mode = "add"
        elif event.key == "control":
            self.mode = "subtract"
        else:
            self.mode = "replace"
            self.selected_bins.clear()
            self._clear_artists(self.selected_artists)

    def _on_motion(self, event):
        """Handle mouse motion event."""
        if event.button != 1:
            return
        if self.drag_start is None or event.inaxes != self.ax:
            return
        if event.xdata is None or event.ydata is None:
            return

        theta = unwrap_theta(self._last_theta, event.xdata)
        self._last_theta = theta

        bins = self._get_bins_from_drag(
            self.drag_start[0], self.drag_start[1],
            theta, event.ydata
        )
        self._update_hover_display(bins)

    def _on_release(self, event):
        """Handle mouse release event."""
        if self.drag_start is None:
            return
        if event.xdata is None or event.ydata is None:
            self.drag_start = None
            self._last_theta = None
            return

        theta = unwrap_theta(self._last_theta, event.xdata)

        bins = self._get_bins_from_drag(
            self.drag_start[0], self.drag_start[1],
            theta, event.ydata
        )

        if self.mode == "replace":
            self.selected_bins.update(bins)
        elif self.mode == "add":
            self.selected_bins.update(bins)
        elif self.mode == "subtract":
            self.selected_bins.difference_update(bins)

        self._update_selected_display()
        self._clear_artists(self.hover_artists)
        self.drag_start = None
        self._last_theta = None

    def _disconnect(self):
        """Disconnect all event handlers."""
        for cid in self._connections:
            self.fig.canvas.mpl_disconnect(cid)
        self.clear_selection()


# Example usage
def create_example_plot():
    """Create an example plot with the polar bin selector."""
    # Data setup
    theta_edges = np.linspace(0, 2 * np.pi, 73)
    r_edges = np.linspace(0, 5, 41)

    theta_centers = 0.5 * (theta_edges[:-1] + theta_edges[1:])
    r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
    theta_vals, r_vals = np.meshgrid(theta_centers, r_centers)
    example = np.sin(2 * theta_vals) * np.exp(-r_vals / 2)

    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, layout='constrained')
    ax.pcolormesh(theta_edges, r_edges, example, shading="auto")
    ax.set_title("Click + Drag = preview | Release = confirm\nShift = add | Ctrl = subtract")
    ax.set_axis_off()

    # Create selector
    selector = SiteSelector(ax, theta_edges, r_edges)

    return fig, ax, selector
