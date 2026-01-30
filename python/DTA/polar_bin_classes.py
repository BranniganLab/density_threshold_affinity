#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 16:22:10 2026.

@author: js2746
"""

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class BinEdge:
    """
    Geometric description of a single visible bin edge in polar coordinates.

    Attributes
    ----------
    r_endpoints : tuple[float, float]
        Radial coordinates of the edge endpoints.
    theta_endpoints : tuple[float, float]
        Angular coordinates (radians) of the edge endpoints.
    """

    r_endpoints: tuple
    theta_endpoints: tuple


class PolarBinGrid:
    """
    Geometry and topology of a polar bin grid.

    This class provides purely computational functionality:
    - Mapping points to bins
    - Determining which bins intersect a polar region
    - Computing which bin edges are externally visible

    It contains no rendering logic and no mutable selection state.
    """

    def __init__(self, theta_edges, r_edges):
        """
        Create a polar bin grid.

        Parameters
        ----------
        theta_edges : array-like
            Angular bin edges in radians.
        r_edges : array-like
            Radial bin edges.
        """
        self.theta_edges = np.asarray(theta_edges)
        self.r_edges = np.asarray(r_edges)
        self.n_t = len(theta_edges) - 1
        self.n_r = len(r_edges) - 1

    def map_coord_to_bin_idx(self, r, theta):
        """
        Determine which bin contains a given polar coordinate.

        Parameters
        ----------
        r : float
            Radial coordinate.
        theta : float
            Angular coordinate in radians.

        Returns
        -------
        tuple[int, int] or None
            The (radial index, angular index) of the bin,
            or None if the point lies outside the grid.
        """
        ti = np.searchsorted(self.theta_edges, theta % (2 * np.pi), side="right") - 1
        ri = np.searchsorted(self.r_edges, r, side="right") - 1

        if 0 <= ri < self.n_r and 0 <= ti < self.n_t:
            return ri, ti
        return None

    def bins_in_region(self, r0, theta0, r1, theta1):
        """
        Return all bins intersecting a dragged polar region.

        The region is defined by two polar coordinates. Angular wraparound
        across 0 / 2π is handled correctly.

        Parameters
        ----------
        r0, theta0 : float
            First corner of the region.
        r1, theta1 : float
            Opposite corner of the region.

        Returns
        -------
        list[tuple[int, int]]
            All bin indices intersecting the region.
        """
        r_min, r_max = sorted((r0, r1))
        bins = []

        for ti in range(self.n_t):
            t_low = self.theta_edges[ti]
            t_high = self.theta_edges[ti + 1]
            if self.bin_in_theta_arc(theta0, theta1, t_low, t_high):
                for ri in range(self.n_r):
                    r_low = self.r_edges[ri]
                    r_high = self.r_edges[ri + 1]
                    if r_high >= r_min and r_low <= r_max:
                        bins.append((ri, ti))
        return bins

    def exposed_edges(self, bins):
        """
        Compute all externally visible edges of a set of bins.

        Internal edges shared by adjacent bins are omitted.

        Parameters
        ----------
        bins : iterable of (int, int)
            Bin indices.

        Returns
        -------
        list[BinEdge]
            Visible boundary edges.
        """
        if not bins:
            return []

        mask = np.zeros((self.n_r, self.n_t), dtype=bool)
        for ri, ti in bins:
            mask[ri, ti] = True

        edges = []
        for ri, ti in zip(*np.where(mask)):
            edges.extend(self._determine_exposed_edges(mask, ri, ti))

        return edges

    def _determine_exposed_edges(self, mask, ri, ti):
        """
        Determine which edges of a bin are exposed.

        Parameters
        ----------
        mask : ndarray
            Boolean array indicating selected bins.
        ri, ti : int
            Bin indices.

        Returns
        -------
        list[BinEdge]
            Exposed edges for the bin.
        """
        edges = []
        n_r, n_t = mask.shape

        if ri == n_r - 1 or not mask[ri + 1, ti]:
            edges.append(self._edge_geometry(ri, ti, "outer"))
        if ri == 0 or not mask[ri - 1, ti]:
            edges.append(self._edge_geometry(ri, ti, "inner"))
        if not mask[ri, (ti - 1) % n_t]:
            edges.append(self._edge_geometry(ri, ti, "left"))
        if not mask[ri, (ti + 1) % n_t]:
            edges.append(self._edge_geometry(ri, ti, "right"))

        return edges

    def _edge_geometry(self, ri, ti, side):
        """
        Construct the geometry for a specific edge of a bin.

        Parameters
        ----------
        ri, ti : int
            Bin indices.
        side : {'outer', 'inner', 'left', 'right'}

        Returns
        -------
        BinEdge
            Edge geometry.
        """
        r0, r1 = self.r_edges[ri], self.r_edges[ri + 1]
        t0 = self.theta_edges[ti]
        t1 = self.theta_edges[(ti + 1) % self.n_t]

        if side == "outer":
            return BinEdge((r1, r1), (t0, t1))
        if side == "inner":
            return BinEdge((r0, r0), (t0, t1))
        if side == "left":
            return BinEdge((r0, r1), (t0, t0))
        if side == "right":
            return BinEdge((r0, r1), (t1, t1))

        raise ValueError(f"Unknown edge type: {side}")

    # pylint: disable=no-else-return
    def bin_in_theta_arc(self, theta_start, theta_end, bin_start, bin_end):
        """Determine if theta arc contains a particular bin.

        Return True if the bin [bin_start, bin_end] intersects the directed angular
        interval from theta_start to theta_end.

        Parameters
        ----------
        theta_start : float
            The starting value for the theta arc.
        theta_end : float
            The ending value for the theta arc.
        bin_start : float
            The starting value for the bin in question.
        bin_end : float
            The ending value for the bin in question.

        Returns
        -------
        Boolean
        """
        TWO_PI = 2 * np.pi

        # Directed arc length
        dtheta = theta_end - theta_start

        # Full circle selects everything
        if abs(dtheta) >= TWO_PI:
            return True

        # Normalize arc start only
        theta_start = theta_start % TWO_PI
        theta_end = theta_start + dtheta

        # Bin edges
        bin_start_n = bin_start % TWO_PI
        if bin_end <= TWO_PI:
            # DO NOT modulo bin_end if it is exactly 2π (or less)
            bin_end_n = bin_end
        else:
            bin_end_n = bin_end % TWO_PI

        if dtheta >= 0:
            # CCW arc
            if theta_end <= TWO_PI:
                # no wrap
                return not (bin_end_n <= theta_start or bin_start_n >= theta_end)
            else:
                # wraps past 2π
                return not (
                    (bin_end_n <= theta_start) and
                    (bin_start_n >= (theta_end - TWO_PI))
                )
        else:
            # CW arc
            if theta_end >= 0:
                # no wrap
                return not (bin_end_n <= theta_end or bin_start_n >= theta_start)
            else:
                # wraps below 0
                return not (
                    (bin_end_n <= (theta_end + TWO_PI)) and
                    (bin_start_n >= theta_start)
                )


class PolarBinRenderer:
    """
    Renderer for drawing polar bin edges on a matplotlib Axes.

    This class converts internal (r, theta) geometry into the
    coordinate order expected by matplotlib polar plots.
    """

    def __init__(self, ax, plot_kwargs=None):
        """
        Create a PolarBinRenderer object and tie it to an Axes instance.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Polar axes used for drawing.
        plot_kwargs : dict
            Dictionary of matplotlib.pyplot keywords for drawing bin edges.
        """
        self.ax = ax
        default_plot_kwargs = {
            "color": "red",
            "lw": 2.0,
            "zorder": 20
        }
        if plot_kwargs is None:
            self.plot_kwargs = default_plot_kwargs
        else:
            self.plot_kwargs = default_plot_kwargs | plot_kwargs

    def draw_edges(self, edges, plot_kwargs=None):
        """
        Draw a collection of bin edges.

        Parameters
        ----------
        edges : iterable of BinEdge
            Edges to draw.
        plot_kwargs : dict
            Keyword arguments passed to ``Axes.plot``.

        Returns
        -------
        list
            Matplotlib artist objects created.
        """
        artists = []
        kwargs = dict(self.plot_kwargs)
        if plot_kwargs is not None:
            kwargs.update(plot_kwargs)

        for edge in edges:
            artists.append(self.ax.plot(edge.theta_endpoints, edge.r_endpoints, **kwargs)[0])
        return artists

    def shade_interior_region(self):
        """Hold space for future method."""


class BinSelectionModel:
    """
    Mutable domain model representing a set of selected bins.

    This class contains no rendering logic and no interaction logic.
    All selection mutations are centralized here to enable future
    undo/redo support.
    """

    def __init__(self):
        """Initialize an empty selection."""
        self._bins = set()

    def snapshot(self):
        """
        Capture the current selection state.

        Returns
        -------
        frozenset
            Immutable snapshot of selected bins.
        """
        return frozenset(self._bins)

    def set(self, bins):
        """
        Replace the current selection.

        Parameters
        ----------
        bins : iterable of (int, int)
            New selection.
        """
        self._bins = set(bins)

    def add(self, bins):
        """
        Add bins to the current selection.

        Parameters
        ----------
        bins : iterable of (int, int)
        """
        self._bins |= set(bins)

    def remove(self, bins):
        """
        Remove bins from the current selection.

        Parameters
        ----------
        bins : iterable of (int, int)
        """
        self._bins -= set(bins)

    def clear(self):
        """Clear the selection."""
        self._bins.clear()

    def bins(self):
        """
        Return the current selection.

        Returns
        -------
        set[tuple[int, int]]
        """
        return set(self._bins)
