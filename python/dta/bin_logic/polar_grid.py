#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Polar bin grid geometry and indexing.

This module defines immutable geometric and topological primitives for
polar binning, including angular and radial bin edges and their indexing
relationships.

The grid:
- provides bin edge definitions and bin index mapping,
- supports geometric queries needed by selection and rendering,
- contains no selection state and no GUI or plotting logic.

It is pure domain code and is shared by GUI and analysis layers.
"""
import itertools
from collections.abc import Iterable
import numpy as np
from dta.bin_logic.utils import Coordinate, BinAddress, BinEdge


class PolarBinGrid:
    """
    Geometry and topology of a polar bin grid.

    This class provides purely computational functionality:
    - Mapping points to bins
    - Determining which bins intersect a polar region
    - Computing which bin edges are externally visible

    It contains no rendering logic and no mutable selection state.
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
        self.d_r = (r_max - r_min) / n_r
        self.d_theta = (2 * np.pi) / n_theta
        self.r_edges = np.linspace(r_min, r_max, n_r + 1)
        self.theta_edges = np.linspace(0.0, 2.0 * np.pi, n_theta + 1)
        self.theta_grid, self.r_grid = np.meshgrid(self.theta_edges, self.r_edges)

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
        theta_idx = (coord[1] % (2 * np.pi)) // self.d_theta
        r_idx = coord[0] // self.d_r

        if 0 <= r_idx < self.n_r and 0 <= theta_idx < self.n_theta:
            return BinAddress(r_idx, theta_idx)
        return None

    def bins_in_region(self,
                       corner1: Coordinate,
                       corner2: Coordinate,
                       span_two_pi: bool = False,
                       ) -> set[BinAddress]:
        """
        Return all bins intersecting a polar region.

        The region is defined by two polar coordinates. Angular wraparound
        across 0 / 2π is handled correctly.

        Parameters
        ----------
        corner1 : Coordinate
            Coordinate of first corner of the region.
        corner2 : Coordinate
            Coordinate of opposite corner of the region.
        span_two_pi : bool, optional
            Whether the region spans two pi or not.

        Returns
        -------
        set[BinAddress]
            All bins intersecting the region.
        """
        start_bin = int(self.map_coord_to_bin_idx(corner1))
        end_bin = int(self.map_coord_to_bin_idx(corner2))

        start_r_index, end_r_index = sorted((start_bin[0], end_bin[0]))
        r_indices = list(range(start_r_index, end_r_index + 1))

        if not span_two_pi:
            start_theta_index, end_theta_index = sorted((start_bin[1], end_bin[1]))
            theta_indices = list(range(start_theta_index, end_theta_index + 1))
        else:
            start_theta_index, end_theta_index = (start_bin[1], end_bin[1])
            if start_theta_index > end_theta_index:
                theta_indices = list(range(0, end_theta_index + 1))
                theta_indices.extend(list(range(start_theta_index, self.n_theta)))
            else:
                theta_indices = list(range(0, start_theta_index + 1))
                theta_indices.extend(list(range(end_theta_index, self.n_theta)))

        bins = list(itertools.product(r_indices, theta_indices))
        return set(bins)

    def exposed_edges(self, bins: Iterable[BinAddress]) -> list[BinEdge]:
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
        if not bins:
            return []

        mask = np.zeros((self.n_r, self.n_theta), dtype=bool)
        for bin_address in bins:
            mask[bin_address] = True

        edges = []
        for bin_address in zip(*np.where(mask)):
            edges.extend(self._determine_exposed_edges(mask, bin_address))

        return edges

    def _determine_exposed_edges(self, mask: np.ndarray, bin_address: BinAddress) -> list[BinEdge]:
        """
        Determine which edges of a bin are exposed.

        Parameters
        ----------
        mask : ndarray
            Boolean array indicating selected bins.
        bin_address : BinAddress
            Bin indices.

        Returns
        -------
        list[BinEdge]
            Exposed edges for the bin.
        """
        edges = []
        n_r, n_t = mask.shape
        ri, ti = bin_address

        if ri == n_r - 1 or not mask[ri + 1, ti]:
            edges.append(self._edge_geometry(bin_address, "outer"))
        if ri == 0 or not mask[ri - 1, ti]:
            edges.append(self._edge_geometry(bin_address, "inner"))
        if not mask[ri, (ti - 1) % n_t]:
            edges.append(self._edge_geometry(bin_address, "left"))
        if not mask[ri, (ti + 1) % n_t]:
            edges.append(self._edge_geometry(bin_address, "right"))

        return edges

    def _edge_geometry(self, bin_address: BinAddress, side: str) -> BinEdge:
        """
        Construct the geometry for a specific edge of a bin.

        Parameters
        ----------
        bin_address : BinAddress
            The r and theta indices of this bin in the lattice.
        side : {'outer', 'inner', 'left', 'right'}

        Returns
        -------
        BinEdge
            Edge geometry.
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

    # pylint: disable=no-else-return
    def bin_in_theta_arc(self,
                         theta_start: float,
                         theta_end: float,
                         bin_start,
                         bin_end,
                         ) -> bool:
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
