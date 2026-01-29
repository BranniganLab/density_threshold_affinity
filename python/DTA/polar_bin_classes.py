#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 16:22:10 2026.

@author: js2746
"""

from dataclasses import dataclass
import numpy as np
from DTA.utils import bin_in_theta_arc


@dataclass(frozen=True)
class BinEdge:
    """
    Geometry for a single exposed bin edge.

    Parameters
    ----------
    r : tuple[float, float]
        Radial coordinates of the edge endpoints.
    theta : tuple[float, float]
        Angular coordinates of the edge endpoints.

    Notes
    -----
    Internally stored as (r, theta) for consistency.
    Converted to (theta, r) only at draw time.
    """
    r: tuple
    theta: tuple


class PolarBinGrid:
    """
    Pure topology and geometry for a polar bin grid.

    This class is responsible for:
    - Mapping points to bins
    - Determining bins intersecting a drag region
    - Computing exposed bin edges

    It contains *no* matplotlib or interaction logic.

    Conventions
    -----------
    - Coordinates: (r, theta)
    - Bin indices: (ri, ti)
    - Masks: mask[ri, ti]
    """

    def __init__(self, theta_edges, r_edges):
        """
        Create a PolarBinGrid.

        Parameters
        ----------
        theta_edges : array-like
            Angular bin edges (radians).
        r_edges : array-like
            Radial bin edges.
        """
        self.theta_edges = np.asarray(theta_edges)
        self.r_edges = np.asarray(r_edges)
        self.n_t = len(theta_edges) - 1
        self.n_r = len(r_edges) - 1

    def bin_at(self, r, theta):
        """
        Return the bin containing a point.

        Parameters
        ----------
        r : float
            Radial coordinate.
        theta : float
            Angular coordinate (radians).

        Returns
        -------
        tuple[int, int] or None
            Bin index (ri, ti), or None if out of bounds.
        """
        ti = np.searchsorted(
            self.theta_edges, theta % (2 * np.pi), side="right"
        ) - 1
        ri = np.searchsorted(self.r_edges, r, side="right") - 1

        if 0 <= ri < self.n_r and 0 <= ti < self.n_t:
            return ri, ti
        return None

    def bins_in_region(self, r_start, theta_start, r_end, theta_end):
        """
        Return all bins intersecting a dragged polar region.

        Angular wraparound (0 / 2π) is handled correctly.

        Parameters
        ----------
        r_start, r_end : float
            Radial drag endpoints.
        theta_start, theta_end : float
            Angular drag endpoints (radians).

        Returns
        -------
        list[tuple[int, int]]
            List of (ri, ti) bin indices.
        """
        r0, r1 = sorted([r_start, r_end])
        bins = []

        for ti in range(self.n_t):
            t_low = self.theta_edges[ti]
            t_high = self.theta_edges[ti + 1]
            if bin_in_theta_arc(theta_start, theta_end, t_low, t_high):
                for ri in range(self.n_r):
                    r_low = self.r_edges[ri]
                    r_high = self.r_edges[ri + 1]
                    if r_high >= r0 and r_low <= r1:
                        bins.append((ri, ti))

        return bins

    def build_mask(self, bins):
        """
        Build a boolean mask from a set of bins.

        Parameters
        ----------
        bins : iterable of (ri, ti)

        Returns
        -------
        ndarray
            Boolean array of shape (n_r, n_t).
        """
        mask = np.zeros((self.n_r, self.n_t), dtype=bool)
        for ri, ti in bins:
            mask[ri, ti] = True
        return mask

    def exposed_edges(self, selected_bins):
        """
        Compute all externally visible edges of a bin set.

        Internal edges shared by adjacent selected bins
        are omitted.

        Parameters
        ----------
        selected_bins : iterable of (ri, ti)

        Returns
        -------
        list[BinEdge]
        """
        if not selected_bins:
            return []

        mask = self.build_mask(selected_bins)
        edges = []

        for ri, ti in zip(*np.where(mask)):
            edges.extend(self._exposed_edges_for_bin(mask, ri, ti))

        return edges

    def _exposed_edges_for_bin(self, mask, ri, ti):
        """
        Determine which edges of a single bin are exposed.

        Parameters
        ----------
        mask : ndarray
            Selection mask.
        ri, ti : int
            Bin indices.

        Returns
        -------
        list[BinEdge]
        """
        edges = []
        n_r, n_t = mask.shape

        if ri == n_r - 1 or not mask[ri + 1, ti]:
            edges.append(self._edge_geometry(ri, ti, "top"))
        if ri == 0 or not mask[ri - 1, ti]:
            edges.append(self._edge_geometry(ri, ti, "bottom"))
        if not mask[ri, (ti - 1) % n_t]:
            edges.append(self._edge_geometry(ri, ti, "left"))
        if not mask[ri, (ti + 1) % n_t]:
            edges.append(self._edge_geometry(ri, ti, "right"))

        return edges

    def _edge_geometry(self, ri, ti, side):
        """
        Return geometry for a single bin edge.

        Parameters
        ----------
        ri, ti : int
            Bin indices.
        side : {'top', 'bottom', 'left', 'right'}

        Returns
        -------
        BinEdge
        """
        t0 = self.theta_edges[ti]
        t1 = self.theta_edges[(ti + 1) % self.n_t]
        r0 = self.r_edges[ri]
        r1 = self.r_edges[ri + 1]

        if side == "top":
            return BinEdge((r1, r1), (t0, t1))
        if side == "bottom":
            return BinEdge((r0, r0), (t0, t1))
        if side == "left":
            return BinEdge((r0, r1), (t0, t0))
        if side == "right":
            return BinEdge((r0, r1), (t1, t1))

        raise ValueError(f"Unknown edge side: {side}")


class PolarBinRenderer:
    """
    Matplotlib renderer for polar bin edges.

    Converts internal (r, theta) geometry into
    matplotlib's (theta, r) plotting convention.
    """

    def __init__(self, ax):
        """
        Create a PolarBinRenderer.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Polar axes to draw on.
        """
        self.ax = ax

    def draw_edges(self, edges, **plot_args):
        """
        Draw a collection of bin edges.

        Parameters
        ----------
        edges : iterable of BinEdge
        **plot_args
            Passed to ``Axes.plot``.

        Returns
        -------
        list of matplotlib.artist.Artist
        """
        artists = []
        for edge in edges:
            artists.append(
                self.ax.plot(edge.theta, edge.r, **plot_args)[0]
            )
        return artists
