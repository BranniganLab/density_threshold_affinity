#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utilities for describing and querying regularly binned polar grids.

This module defines geometry objects used to represent polar coordinate
lattices composed of discrete radial and angular bins. The primary object,
:class:`PolarBinGrid`, models a polar domain partitioned into regularly spaced
radial and angular intervals and provides geometric operations on that grid.

The grid itself stores only lattice geometry. It does not store selections,
perform rendering, or depend on GUI backends. Geometry utilities provided by
this module support analysis workflows, plotting code, and interactive tools
that operate on polar bin structures.

Contents
--------
GridDim
    Immutable description of a single discretized dimension. Stores bounds,
    bin count, derived bin width, and bin edge coordinates.

PolarBinGrid
    Polar coordinate lattice supporting coordinate-to-bin mapping, region
    queries, exposed-edge determination, and bin geometry calculations.

Examples
--------
Create a polar grid with radial bounds from 0 to 5 and 72 angular bins::

    grid = PolarBinGrid(
        r_min=0.0,
        r_max=5.0,
        n_r=40,
        n_theta=72,
    )

Determine which lattice bin contains a coordinate::

    address = grid.map_coord_to_bin_idx((2.5, np.pi))

Find all bins touched by a drag-selection region::

    bins = grid.get_bins_in_region(
        corner1=(1.0, 0.1),
        corner2=(3.0, 0.5),
    )
"""
from __future__ import annotations
import itertools
from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Literal
import numpy as np
from dta.bin_logic.utils import Coordinate, BinAddress, BinEdge

BinSide = Literal["outer", "inner", "left", "right"]


@dataclass(frozen=True)
class GridDim:
    """
    Immutable description of a discretized coordinate dimension.

    A ``GridDim`` defines one regularly spaced dimension of a lattice by
    storing lower and upper bounds together with the number of bins spanning
    that interval. Derived geometric quantities including bin width and bin
    edge coordinates are calculated automatically during initialization.

    ``GridDim`` is intended to separate dimension-specific geometry from
    higher-level grid logic. For example, a polar grid may contain one
    ``GridDim`` describing radial spacing and another describing angular
    spacing.

    Attributes
    ----------
    lower_bound : float
        Lower boundary of the dimension.

    upper_bound : float
        Upper boundary of the dimension.

    n_bins : int
        Number of regularly spaced bins.

    bin_width : float
        Width of each bin, calculated as
        ``(upper_bound - lower_bound) / n_bins``.

    bin_edges : ndarray
        Array containing bin boundary coordinates. Length is
        ``n_bins + 1``.
    """

    lower_bound: float
    upper_bound: float
    n_bins: int
    bin_width: float = field(init=False)
    bin_edges: np.ndarray = field(init=False)

    def __post_init__(self):
        """Calculate bin width and edges from input variables."""
        if not isinstance(self.n_bins, int):
            raise TypeError(f"Expected int but got {type(self.n_bins)} ({self.n_bins}) for n_bins")
        bin_width = (self.upper_bound - self.lower_bound) / self.n_bins
        object.__setattr__(self, 'bin_width', bin_width)
        edges = np.linspace(self.lower_bound, self.upper_bound, self.n_bins + 1)
        object.__setattr__(self, 'bin_edges', edges)

    def __eq__(self, other: GridDim) -> bool | NotImplemented:
        """Compare two GridDims for equality."""
        if not isinstance(other, GridDim):
            return NotImplemented
        return np.array_equal(self.bin_edges, other.bin_edges)


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
    r : GridDim
        The radial dimension of the polar lattice.
    theta : GridDim
        The angular dimension of the polar lattice.
    r_grid : ndarray
        Meshgrid array of radial edge coordinates suitable for plotting with
        polar grid data.
    theta_grid : ndarray
        Meshgrid array of angular edge coordinates suitable for plotting with
        polar grid data.
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
        if r_max <= r_min:
            raise ValueError("r_max must be greater than r_min.")
        if n_r <= 0:
            raise ValueError("n_r must be positive.")
        if n_theta <= 0:
            raise ValueError("n_theta must be positive.")
        self.r = GridDim(r_min, r_max, n_r)
        self.theta = GridDim(0, 2 * np.pi, n_theta)
        self.theta_grid, self.r_grid = np.meshgrid(self.theta.bin_edges, self.r.bin_edges)

    def __eq__(self, other: PolarBinGrid) -> bool | NotImplemented:
        """Compare PolarBinGrid objects for equality."""
        if not isinstance(other, PolarBinGrid):
            return NotImplemented
        if ((self.r == other.r) and (self.theta == other.theta)):
            return True
        return False

    @property
    def bin_areas(self) -> np.ndarray:
        """Calculate area of every bin in lattice."""
        areas = np.zeros((self.r.n_bins, self.theta.n_bins))
        for radial_ring in range(self.r.n_bins):
            bin_address = BinAddress(radial_ring, 0)
            areas[radial_ring, :] = self.calc_bin_area(bin_address)
        return areas

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
        theta_idx = int((coord[1] % (2.0 * np.pi)) // self.theta.bin_width)
        r_idx = int((coord[0] - self.r.lower_bound) // self.r.bin_width)

        if 0 <= r_idx < self.r.n_bins and 0 <= theta_idx < self.theta.n_bins:
            return BinAddress(r_idx, theta_idx)
        return None

    def get_bins_in_region(self,
                           corner1: Coordinate,
                           corner2: Coordinate,
                           crosses_theta_boundary: bool = False,
                           ) -> set[BinAddress]:
        """
        Return all bins intersecting a polar region.

        The region is defined by two polar coordinates. Angular wraparound
        across 0 / 2π is controlled by crosses_theta_boundary.

        Parameters
        ----------
        corner1 : Coordinate
            Coordinate of first corner of the region.
        corner2 : Coordinate
            Coordinate of opposite corner of the region.
        crosses_theta_boundary : bool, optional
            Whether the directed angular interval crosses the periodic 0 / 2π
            boundary. Default is False.

        Returns
        -------
        set[BinAddress]
            All bins intersecting the region.
        """
        corner1_bin = self.map_coord_to_bin_idx(corner1)
        corner2_bin = self.map_coord_to_bin_idx(corner2)

        if not (corner1_bin and corner2_bin):
            raise ValueError("One or both corners were outside of the grid domain.")

        # Drag order doesn't matter. Sort the r values from low to high.
        r_index1, r_index2 = sorted((corner1_bin[0], corner2_bin[0]))
        r_indices = list(range(r_index1, r_index2 + 1))

        if not crosses_theta_boundary:
            # Treat this as a single rectangular region. Drag order doesn't matter.
            theta_index1, theta_index2 = sorted((corner1_bin[1], corner2_bin[1]))
            theta_indices = list(range(theta_index1, theta_index2 + 1))
        else:
            # Treat this as two rectangular regions. Drag order matters.
            theta_index1, theta_index2 = corner1_bin[1], corner2_bin[1]
            if theta_index1 <= theta_index2:
                # Select from 0 to theta_index1 and from theta_index2 to 2pi
                theta_indices = list(range(0, theta_index1 + 1))
                theta_indices.extend(range(theta_index2, self.theta.n_bins))
            else:
                # Select from 0 to theta_index2 and from theta_index1 to 2pi
                theta_indices = list(range(theta_index1, self.theta.n_bins))
                theta_indices.extend(range(0, theta_index2 + 1))

        # Calculate Cartesian product to produce all ordered pairs
        bin_pairs = set(itertools.product(r_indices, theta_indices))

        # Convert to set of BinAddress objects
        bins = {BinAddress(r_idx, theta_idx) for r_idx, theta_idx in bin_pairs}
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
        bins = set(bins)
        if not bins:
            return []

        mask = np.zeros((self.r.n_bins, self.theta.n_bins), dtype=bool)
        for bin_address in bins:
            mask[bin_address] = True

        edges = []
        for bin_address in bins:
            edges.extend(self._determine_exposed_bin_edges(mask, bin_address))

        return edges

    def _determine_exposed_bin_edges(self,
                                     mask: np.ndarray,
                                     bin_address: BinAddress
                                     ) -> list[BinEdge]:
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
            # The bin to the "left" is not selected
            edges.append(self._determine_bin_edge(bin_address, "left"))
        if not mask[ri, (ti + 1) % n_t]:
            # The bin to the "right" is not selected
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
        r0 = self.r.bin_edges[bin_address[0]]
        r1 = self.r.bin_edges[bin_address[0] + 1]
        t0 = self.theta.bin_edges[bin_address[1]]
        t1 = self.theta.bin_edges[(bin_address[1] + 1) % self.theta.n_bins]

        if side == "outer":
            return BinEdge(Coordinate(r1, t0), Coordinate(r1, t1))
        if side == "inner":
            return BinEdge(Coordinate(r0, t0), Coordinate(r0, t1))
        if side == "left":
            return BinEdge(Coordinate(r0, t0), Coordinate(r1, t0))
        if side == "right":
            return BinEdge(Coordinate(r0, t1), Coordinate(r1, t1))

        raise ValueError(f"Unknown edge type: {side}")

    def calc_bin_area(self, bin_address: BinAddress) -> float:
        r"""
        Calculate bin area.

        Formula for bin area in polar lattice is
        $r_i * \delta r * \delta \theta$,
        where $r_i$ is the radial distance to the midpoint of the bin.

        Parameters
        ----------
        bin_address : BinAddress
            The r and theta indices of this bin in the lattice.

        Returns
        -------
        float
            The area of the bin.
        """
        bin_lower_r_bound = bin_address[0] * self.r.bin_width + self.r.lower_bound
        r_i = bin_lower_r_bound + (self.r.bin_width * 0.5)
        return r_i * self.r.bin_width * self.theta.bin_width
