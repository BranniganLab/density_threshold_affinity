#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for PolarBinGrid.
"""

import numpy as np
import pytest

from dta.bin_logic import PolarBinGrid
from dta.bin_logic.utils import BinAddress


def test_init_sets_expected_grid_attributes():
    """
    Test PolarBinGrid initialization.

    Test that it:
        1) stores radial and angular bin counts
        2) stores radial and angular bin widths
        3) creates radial and angular edge arrays with expected values
        4) creates meshgrid arrays with expected shapes
    """
    grid = PolarBinGrid(1.0, 3.0, 2, 4)

    # stores bin counts
    assert grid.n_r == 2
    assert grid.n_theta == 4

    # stores bin widths
    assert grid.d_r == 1.0
    assert grid.d_theta == 0.5 * np.pi

    # creates expected bin edge arrays
    np.testing.assert_allclose(grid.r_edges, np.array([1.0, 2.0, 3.0]))
    np.testing.assert_allclose(
        grid.theta_edges,
        np.array([0.0, 0.5 * np.pi, np.pi, 1.5 * np.pi, 2.0 * np.pi]),
    )

    # creates meshgrid arrays with one entry for each r/theta edge pair
    assert grid.r_grid.shape == (3, 5)
    assert grid.theta_grid.shape == (3, 5)


def test_init_rejects_invalid_grid_parameters():
    """
    Test PolarBinGrid input validation.

    Test that it:
        1) rejects radial bounds with no positive width
        2) rejects zero or negative radial bin counts
        3) rejects zero or negative angular bin counts
    """
    # rejects radial bounds with no positive width
    with pytest.raises(ValueError, match="r_max must be greater than r_min"):
        PolarBinGrid(1.0, 1.0, 2, 4)

    # rejects zero or negative radial bin counts
    with pytest.raises(ValueError, match="n_r must be positive"):
        PolarBinGrid(0.0, 1.0, 0, 4)

    # rejects zero or negative angular bin counts
    with pytest.raises(ValueError, match="n_theta must be positive"):
        PolarBinGrid(0.0, 1.0, 2, 0)


def test_map_coord_to_bin_idx():
    """
    Test map_coord_to_bin_idx.

    Test that it:
        1) maps to the correct bin
        2) handles wrapping around 2pi correctly
        3) ignores values outside of lattice boundaries
        4) handles coords on bin boundary correctly (default to higher bin index)
    """
    grid = PolarBinGrid(0, 1, 2, 4)

    r = 0.25
    theta = 0.25 * np.pi

    # maps to correct bin
    assert grid.map_coord_to_bin_idx((r, theta)) == BinAddress(0, 0)

    # correctly wraps around 2pi
    assert grid.map_coord_to_bin_idx((r, theta + 2.0 * np.pi)) == BinAddress(0, 0)

    # outside lattice boundary is None
    assert grid.map_coord_to_bin_idx((2.0, theta)) is None

    # defaults to higher bin index when coord on bin boundary
    assert grid.map_coord_to_bin_idx((0.5, np.pi)) == BinAddress(1, 2)

    # exactly on outer edge is None because it defaults to higher bin index
    assert grid.map_coord_to_bin_idx((1.0, theta)) is None


def test_map_coord_to_bin_idx_with_nonzero_r_min():
    """
    Test map_coord_to_bin_idx with a nonzero radial minimum.

    Test that it:
        1) uses r_min when calculating radial bin index
        2) maps the first radial interval correctly
        3) maps the second radial interval correctly
        4) rejects coordinates below r_min
    """
    grid = PolarBinGrid(1.0, 3.0, 2, 4)

    theta = 0.25 * np.pi

    # uses r_min when calculating the first radial bin
    assert grid.map_coord_to_bin_idx((1.25, theta)) == BinAddress(0, 0)

    # maps the upper radial interval to the second radial bin
    assert grid.map_coord_to_bin_idx((2.25, theta)) == BinAddress(1, 0)

    # defaults to higher bin index on radial boundary
    assert grid.map_coord_to_bin_idx((2.0, theta)) == BinAddress(1, 0)

    # below r_min is outside the grid
    assert grid.map_coord_to_bin_idx((0.75, theta)) is None


def test_get_bins_in_region_non_wrapping_region():
    """
    Test get_bins_in_region without angular wraparound.

    Test that it:
        1) returns all bins between the two radial indices
        2) returns all bins between the two angular indices
        3) returns BinAddress objects rather than plain tuples
    """
    grid = PolarBinGrid(0.0, 1.0, 4, 8)

    bins = grid.get_bins_in_region((0.1, 0.1), (0.6, 0.6), False)

    expected = {
        BinAddress(0, 0),
        BinAddress(1, 0),
        BinAddress(2, 0),
    }

    # returns all bins in the rectangular polar region
    assert bins == expected

    # returns BinAddress objects
    assert all(isinstance(bin_address, BinAddress) for bin_address in bins)


def test_get_bins_in_region_reversed_corners():
    """
    Test get_bins_in_region with corners supplied in reverse order.

    Test that it:
        1) handles decreasing radial coordinates
        2) handles decreasing angular coordinates
        3) returns the same bin set as the equivalent forward-order region
    """
    grid = PolarBinGrid(0.0, 1.0, 4, 8)

    forward_bins = grid.get_bins_in_region((0.1, 0.1), (0.6, 0.6), False)
    reversed_bins = grid.get_bins_in_region((0.6, 0.6), (0.1, 0.1), False)

    # corner order does not affect a non-wrapping region
    assert reversed_bins == forward_bins


def test_get_bins_in_region_wraparound_includes_boundary_bins():
    """
    Test get_bins_in_region with angular wraparound.

    Test that it:
        1) includes bins near theta equals zero
        2) includes bins near theta equals 2pi
        3) excludes bins in the unselected middle angular interval
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 8)

    bins = grid.get_bins_in_region(
        (0.1, 1.9 * np.pi),
        (0.9, 0.1 * np.pi),
        True,
    )

    # includes bins near theta equals zero
    assert grid.map_coord_to_bin_idx((0.5, 0.01)) in bins

    # includes bins near theta equals 2pi
    assert grid.map_coord_to_bin_idx((0.5, 1.99 * np.pi)) in bins

    # excludes bins in the middle angular interval
    assert grid.map_coord_to_bin_idx((0.5, np.pi)) not in bins


def test_get_bins_in_region_wraparound_matching_theta_bins_selects_all_theta():
    """
    Test get_bins_in_region when a wrapped region starts and ends in same theta bin.

    Test that it:
        1) treats matching theta bins as a full angular span
        2) includes every angular bin
        3) includes the requested radial range
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 8)

    bins = grid.get_bins_in_region((0.1, 0.1), (0.9, 0.2), True)

    expected = {
        BinAddress(r_idx, theta_idx)
        for r_idx in range(grid.n_r)
        for theta_idx in range(grid.n_theta)
    }

    # selecting across 2pi with matching theta bins selects all angular bins
    assert bins == expected


def test_get_bins_in_region_rejects_out_of_bounds_corners():
    """
    Test get_bins_in_region with corners outside the grid.

    Test that it:
        1) rejects a first corner outside the radial domain
        2) rejects a second corner outside the radial domain
        3) raises a clear ValueError rather than failing from None indexing
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 4)

    # rejects a first corner outside the grid
    with pytest.raises(ValueError, match="outside of the grid domain"):
        grid.get_bins_in_region((2.0, 0.1), (0.5, 0.2), False)

    # rejects a second corner outside the grid
    with pytest.raises(ValueError, match="outside of the grid domain"):
        grid.get_bins_in_region((0.5, 0.1), (2.0, 0.2), False)


def test_list_all_exposed_edges_empty_input_has_no_edges():
    """
    Test list_all_exposed_edges with no bins.

    Test that it:
        1) accepts an empty iterable
        2) returns an empty list
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 4)

    edges = grid.list_all_exposed_edges([])

    # empty bin collections have no exposed edges
    assert edges == []


def test_list_all_exposed_edges_single_bin_has_four_edges():
    """
    Test list_all_exposed_edges with a single bin.

    Test that it:
        1) returns one edge for each side of the bin
        2) includes inner and outer radial edges
        3) includes left and right angular edges
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 4)

    edges = grid.list_all_exposed_edges({BinAddress(0, 0)})

    # one isolated bin has four exposed edges
    assert len(edges) == 4


def test_list_all_exposed_edges_two_theta_adjacent_bins_share_internal_edge():
    """
    Test list_all_exposed_edges with two angularly adjacent bins.

    Test that it:
        1) omits the internal edge shared by adjacent theta bins
        2) preserves the six externally visible edges
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 4)

    edges = grid.list_all_exposed_edges({BinAddress(0, 0), BinAddress(0, 1)})

    # two adjacent bins share one internal edge, removing two exposed edges
    assert len(edges) == 6


def test_list_all_exposed_edges_two_radially_adjacent_bins_share_internal_edge():
    """
    Test list_all_exposed_edges with two radially adjacent bins.

    Test that it:
        1) omits the internal edge shared by adjacent radial bins
        2) preserves the six externally visible edges
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 4)

    edges = grid.list_all_exposed_edges({BinAddress(0, 0), BinAddress(1, 0)})

    # two radially adjacent bins share one internal edge, removing two exposed edges
    assert len(edges) == 6


def test_list_all_exposed_edges_theta_wraparound_bins_share_internal_edge():
    """
    Test list_all_exposed_edges with bins adjacent across the angular boundary.

    Test that it:
        1) treats theta index zero and the final theta index as neighbors
        2) omits their shared internal edge
        3) preserves the six externally visible edges
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 4)

    edges = grid.list_all_exposed_edges({BinAddress(0, 0), BinAddress(0, 3)})

    # bins adjacent across theta equals zero share one internal edge
    assert len(edges) == 6


def test_determine_bin_edge_returns_expected_coordinates():
    """
    Test _determine_bin_edge.

    Test that it:
        1) returns the expected outer edge
        2) returns the expected inner edge
        3) returns the expected left edge
        4) returns the expected right edge
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 4)
    bin_address = BinAddress(0, 0)

    # returns expected outer radial edge
    assert grid._determine_bin_edge(bin_address, "outer") == (
        (0.5, 0.0),
        (0.5, 0.5 * np.pi),
    )

    # returns expected inner radial edge
    assert grid._determine_bin_edge(bin_address, "inner") == (
        (0.0, 0.0),
        (0.0, 0.5 * np.pi),
    )

    # returns expected left angular edge
    assert grid._determine_bin_edge(bin_address, "left") == (
        (0.0, 0.0),
        (0.5, 0.0),
    )

    # returns expected right angular edge
    assert grid._determine_bin_edge(bin_address, "right") == (
        (0.0, 0.5 * np.pi),
        (0.5, 0.5 * np.pi),
    )


def test_determine_bin_edge_rejects_unknown_side():
    """
    Test _determine_bin_edge input validation.

    Test that it:
        1) rejects an unknown edge side
        2) raises a clear ValueError
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 4)

    # rejects invalid edge side names
    with pytest.raises(ValueError, match="Unknown edge type"):
        grid._determine_bin_edge(BinAddress(0, 0), "diagonal")
