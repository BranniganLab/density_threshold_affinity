#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utilities for describing and querying a regularly binned polar grid.

This module defines :class:`PolarBinGrid`, a small geometry object that
represents a polar coordinate grid divided into radial and angular bins. The
grid covers radii from ``r_min`` to ``r_max`` and angles from ``0`` to
``2*pi``. Each bin is addressed by a ``BinAddress`` containing a radial index
and an angular index.

The module provides operations for converting polar coordinates to bin
indices, finding all bins touched by a rectangular polar region, and computing
the exposed boundary edges of a set of bins. These operations are useful for
analysis code, plotting code, and interactive selection tools, but the module
itself does not store selections, draw figures, or depend on a GUI backend.
"""
import itertools
from collections.abc import Iterable
from typing import Literal
import numpy as np
from dta.bin_logic.utils import Coordinate, BinAddress, BinEdge

BinSide = Literal["outer", "inner", "left", "right"]


class PolarBinGrid:
    """
    Regular polar grid with radial and angular bin indexing.

    A ``PolarBinGrid`` divides a polar coordinate domain into ``n_r`` radial
    bins and ``n_theta`` angular bins. Radial bins span the interval from
    ``r_min`` to ``r_max``. Angular bins span the full circle from ``0`` to
    ``2*pi`` and wrap periodically at the angular boundary.

    The class is responsible for geometric and topological queries on that
    grid. It can map an ``(r, theta)`` coordinate to the bin containing it,
    enumerate the bins intersected by a polar region, and identify the outer
    boundary edges of an arbitrary collection of bins. It does not track which
    bins are selected and does not perform any rendering.

    Attributes
    ----------
    n_r : int
        Number of radial bins in the grid.
    n_theta : int
        Number of angular bins in the grid.
    r_edges : ndarray
        Radial bin edge coordinates with shape ``(n_r + 1,)``.
    theta_edges : ndarray
        Angular bin edge coordinates with shape ``(n_theta + 1,)`` spanning
        ``0`` to ``2*pi``.
    r_grid : ndarray
        Meshgrid array of radial edge coordinates suitable for plotting.
    theta_grid : ndarray
        Meshgrid array of angular edge coordinates suitable for plotting.

    Properties
    ----------
    r_min, r_max, d_r, d_theta.
    """

    def __init__(self,
                 r_min: float,
                 r_max: float,
                 n_r: int,
                 n_theta: int,
                 ) -> None:
        """
        Create a polar bin grid.

        Parameters
        ----------
        r_min : float
            The minimum r value in the grid. The beginning of the first bin.
        r_max : float
            The maximum r value in the grid. The end of the last bin.
        n_r : int
            Number of radial (r) bins.
        n_theta : int
            Number of angular (theta) bins.
        """
        self.n_r = n_r
        self.n_theta = n_theta
        self.r_edges = np.linspace(r_min, r_max, n_r + 1)
        self.theta_edges = np.linspace(0.0, 2.0 * np.pi, n_theta + 1)
        self.theta_grid, self.r_grid = np.meshgrid(self.theta_edges, self.r_edges)

    @property
    def r_min(self) -> float:
        """The minimum radial value (beginning of first bin)."""
        return self.r_edges[0]

    @property
    def r_max(self) -> float:
        """The maximum radial value (end of last bin)."""
        return self.r_edges[-1]

    @property
    def d_r(self) -> float:
        """The radial step size dr."""
        return (self.r_max - self.r_min) / self.n_r

    @property
    def d_theta(self) -> float:
        """The angular step size dtheta."""
        return (2 * np.pi) / self.n_theta

    def map_coord_to_bin_idx(self, coord: Coordinate) -> BinAddress | None:
        """
        Determine which bin contains a given polar coordinate.

        Parameters
        ----------
        coord : Coordinate
            The r, theta coordinate pair.

        Returns
        -------
        BinAddress or None
            The (radial index, angular index) of the bin,
            or None if the point lies outside the grid.
        """
        theta_idx = int((coord[1] % (2 * np.pi)) // self.d_theta)
        r_idx = int(coord[0] // self.d_r)

        if 0 <= r_idx < self.n_r and 0 <= theta_idx < self.n_theta:
            return BinAddress(r_idx, theta_idx)
        return None

    def get_bins_in_region(self,
                           corner1: Coordinate,
                           corner2: Coordinate,
                           span_two_pi: bool = False,
                           ) -> set[BinAddress]:
        """
        Return all bins intersecting a polar region.

        The region is defined by two polar coordinates. Angular wraparound
        across 0 / 2π is handled depending on the span_two_pi flag.

        Parameters
        ----------
        corner1 : Coordinate
            Coordinate of first corner of the region.
        corner2 : Coordinate
            Coordinate of opposite corner of the region.
        span_two_pi : bool, optional
            Whether the region spans 2*pi or not.

        Returns
        -------
        set[BinAddress]
            All bins intersecting the region.
        """
        corner1_bin = self.map_coord_to_bin_idx(corner1)
        corner2_bin = self.map_coord_to_bin_idx(corner2)

        r_index1, r_index2 = sorted((corner1_bin[0], corner2_bin[0]))
        r_indices = list(range(r_index1, r_index2 + 1))

        if span_two_pi is False:
            # Treat this as a singular rectangular region
            theta_index1, theta_index2 = sorted((corner1_bin[1], corner2_bin[1]))
            theta_indices = list(range(theta_index1, theta_index2 + 1))
        else:
            # Treat this as two rectangular regions
            theta_index1, theta_index2 = (corner1_bin[1], corner2_bin[1])
            if theta_index1 == theta_index2:
                # Select all theta values
                theta_indices = list(range(0, self.n_theta))
            elif theta_index1 < theta_index2:
                # Select from 0 to theta1 and from theta2 to ntheta
                theta_indices = list(range(0, theta_index1 + 1))
                theta_indices.extend(list(range(theta_index2, self.n_theta)))
            else:
                # Select from 0 to theta2 and from theta1 to ntheta
                theta_indices = list(range(0, theta_index2 + 1))
                theta_indices.extend(list(range(theta_index1, self.n_theta)))

        # Calculate Cartesian product to produce all ordered pairs
        bins = set(itertools.product(r_indices, theta_indices))
        return bins

    def list_all_exposed_edges(self, bins: Iterable[BinAddress]) -> list[BinEdge]:
        """
        Compute all externally visible edges of a set of bins.

        Internal edges shared by adjacent bins are omitted.

        Parameters
        ----------
        bins : iterable of BinAddress(es)
            Bin indices.

        Returns
        -------
        list[BinEdge]
            Visible boundary edges.
        """
        if len(bins) == 0:
            return []

        mask = np.zeros((self.n_r, self.n_theta), dtype=bool)
        for bin_address in bins:
            mask[bin_address] = True

        edges = []
        for bin_address in bins:
            edges.extend(self._determine_exposed_bin_edges(mask, bin_address))

        return edges

    def _determine_exposed_bin_edges(self, mask: np.ndarray, bin_address: BinAddress) -> list[BinEdge]:
        """
        Determine which edges of a bin are exposed.

        Parameters
        ----------
        mask : ndarray
            Boolean array indicating all selected bins.
        bin_address : BinAddress
            The r and theta index of the bin whose edges you wish to check.

        Returns
        -------
        list[BinEdge]
            Exposed edges for the bin.
        """
        edges = []
        n_r, n_t = mask.shape
        ri, ti = bin_address

        if ri == n_r - 1 or not mask[ri + 1, ti]:
            # Bin is in outermost radial shell
            # OR
            # the bin in the next radial shell is not selected
            edges.append(self._determine_bin_edge(bin_address, "outer"))
        if ri == 0 or not mask[ri - 1, ti]:
            # Bin is in innermost radial shell
            # OR
            # the bin in the previous radial shell is not selected
            edges.append(self._determine_bin_edge(bin_address, "inner"))
        if not mask[ri, (ti - 1) % n_t]:
            # The bin to the "left" is empty
            edges.append(self._determine_bin_edge(bin_address, "left"))
        if not mask[ri, (ti + 1) % n_t]:
            # The bin to the "right" is empty
            edges.append(self._determine_bin_edge(bin_address, "right"))

        return edges

    def _determine_bin_edge(self,
                            bin_address: BinAddress,
                            side: BinSide,
                            ) -> BinEdge:
        """
        Compute the coordinates that define a BinEdge from an address and side.

        Parameters
        ----------
        bin_address : BinAddress
            The r and theta indices of this bin in the lattice.
        side : {'outer', 'inner', 'left', 'right'}

        Returns
        -------
        BinEdge
            The two coordinates that define a bin's edge on the polar lattice.
        """
        r0 = self.r_edges[bin_address[0]]
        r1 = self.r_edges[bin_address[0] + 1]
        t0 = self.theta_edges[bin_address[1]]
        t1 = self.theta_edges[(bin_address[1] + 1) % self.n_theta]

        if side == "outer":
            return BinEdge(Coordinate(r1, t0), Coordinate(r1, t1))
        if side == "inner":
            return BinEdge(Coordinate(r0, t0), Coordinate(r0, t1))
        if side == "left":
            return BinEdge(Coordinate(r0, t0), Coordinate(r1, t0))
        if side == "right":
            return BinEdge(Coordinate(r0, t1), Coordinate(r1, t1))

        raise ValueError(f"Unknown edge type: {side}")
