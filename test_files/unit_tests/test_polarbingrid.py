#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for GridDim and PolarBinGrid.
"""

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from dta.bin_logic.polar_grid import GridDim, PolarBinGrid
from dta.bin_logic.utils import BinAddress


# =============================================================================
# GridDim tests
# =============================================================================


def test_griddim_calculates_expected_attributes():
    """
    Test GridDim initialization.

    Test that it:
        1) stores lower and upper bounds
        2) stores the number of bins
        3) calculates bin width
        4) creates bin edge array with expected values
    """
    dim = GridDim(1.0, 3.0, 2)

    assert dim.lower_bound == 1.0
    assert dim.upper_bound == 3.0
    assert dim.n_bins == 2

    np.testing.assert_allclose(dim.bin_width, 1.0)
    np.testing.assert_allclose(
        dim.bin_edges,
        np.array([1.0, 2.0, 3.0]),
    )


def test_griddim_calculates_expected_angular_dimension():
    """
    Test GridDim initialization for an angular dimension.

    Test that it:
        1) accepts zero as the lower bound
        2) accepts 2pi as the upper bound
        3) calculates angular bin width
        4) creates angular bin edge array with expected values
    """
    dim = GridDim(0.0, 2.0 * np.pi, 4)

    assert dim.lower_bound == 0.0
    np.testing.assert_allclose(dim.upper_bound, 2.0 * np.pi)
    assert dim.n_bins == 4

    np.testing.assert_allclose(dim.bin_width, 0.5 * np.pi)
    np.testing.assert_allclose(
        dim.bin_edges,
        np.array([
            0.0,
            0.5 * np.pi,
            np.pi,
            1.5 * np.pi,
            2.0 * np.pi,
        ]),
    )


def test_griddim_is_frozen():
    """
    Test GridDim immutability.

    Test that it:
        1) prevents reassignment of input attributes
        2) prevents reassignment of calculated attributes
    """
    dim = GridDim(1.0, 3.0, 2)

    with pytest.raises(FrozenInstanceError):
        dim.lower_bound = 0.0

    with pytest.raises(FrozenInstanceError):
        dim.bin_width = 10.0


# =============================================================================
# PolarBinGrid initialization tests
# =============================================================================


def test_polarbingrid_init_sets_expected_dimension_attributes():
    """
    Test PolarBinGrid initialization.

    Test that it:
        1) stores radial GridDim attributes
        2) stores angular GridDim attributes
        3) creates radial dimension with expected values
        4) creates angular dimension with expected values
    """
    grid = PolarBinGrid(1.0, 3.0, 2, 4)

    assert isinstance(grid.r, GridDim)
    assert isinstance(grid.theta, GridDim)

    assert grid.r.lower_bound == 1.0
    assert grid.r.upper_bound == 3.0
    assert grid.r.n_bins == 2

    np.testing.assert_allclose(grid.r.bin_width, 1.0)
    np.testing.assert_allclose(
        grid.r.bin_edges,
        np.array([1.0, 2.0, 3.0]),
    )

    assert grid.theta.lower_bound == 0.0
    np.testing.assert_allclose(
        grid.theta.upper_bound,
        2.0 * np.pi,
    )
    assert grid.theta.n_bins == 4

    np.testing.assert_allclose(
        grid.theta.bin_width,
        0.5 * np.pi,
    )


def test_polarbingrid_init_sets_expected_meshgrid_attributes():
    """
    Test PolarBinGrid meshgrid initialization.

    Test that it:
        1) creates radial and angular meshgrid arrays
        2) creates meshgrid arrays with expected shapes
        3) uses radial bin edges in radial meshgrid
        4) uses angular bin edges in angular meshgrid
    """
    grid = PolarBinGrid(1.0, 3.0, 2, 4)

    assert grid.r_grid.shape == (3, 5)
    assert grid.theta_grid.shape == (3, 5)

    np.testing.assert_allclose(
        grid.r_grid[:, 0],
        grid.r.bin_edges,
    )

    np.testing.assert_allclose(
        grid.theta_grid[0, :],
        grid.theta.bin_edges,
    )


def test_polarbingrid_init_rejects_invalid_grid_parameters():
    """
    Test PolarBinGrid input validation.

    Test that it:
        1) rejects radial bounds with no positive width
        2) rejects zero or negative radial bin counts
        3) rejects zero or negative angular bin counts
    """
    with pytest.raises(
        ValueError,
        match="r_max must be greater than r_min",
    ):
        PolarBinGrid(1.0, 1.0, 2, 4)

    with pytest.raises(
        ValueError,
        match="n_r must be positive",
    ):
        PolarBinGrid(0.0, 1.0, 0, 4)

    with pytest.raises(
        ValueError,
        match="n_theta must be positive",
    ):
        PolarBinGrid(0.0, 1.0, 2, 0)


# =============================================================================
# Coordinate mapping tests
# =============================================================================


def test_map_coord_to_bin_idx():
    """
    Test map_coord_to_bin_idx.

    Test that it:
        1) maps to the correct bin
        2) handles wrapping around 2pi correctly
        3) ignores values outside lattice boundaries
        4) handles coords on bin boundary correctly
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 4)

    r = 0.25
    theta = 0.25 * np.pi

    assert grid.map_coord_to_bin_idx(
        (r, theta),
    ) == BinAddress(0, 0)

    assert grid.map_coord_to_bin_idx(
        (r, theta + 2.0 * np.pi),
    ) == BinAddress(0, 0)

    assert grid.map_coord_to_bin_idx(
        (2.0, theta),
    ) is None

    assert grid.map_coord_to_bin_idx(
        (0.5, np.pi),
    ) == BinAddress(1, 2)

    assert grid.map_coord_to_bin_idx(
        (1.0, theta),
    ) is None


def test_map_coord_to_bin_idx_with_nonzero_r_min():
    """
    Test map_coord_to_bin_idx with nonzero radial minimum.

    Test that it:
        1) uses r_min in radial index calculation
        2) maps first radial interval correctly
        3) maps second radial interval correctly
        4) rejects coordinates below r_min
    """
    grid = PolarBinGrid(1.0, 3.0, 2, 4)

    theta = 0.25 * np.pi

    assert grid.map_coord_to_bin_idx(
        (1.25, theta),
    ) == BinAddress(0, 0)

    assert grid.map_coord_to_bin_idx(
        (2.25, theta),
    ) == BinAddress(1, 0)

    assert grid.map_coord_to_bin_idx(
        (2.0, theta),
    ) == BinAddress(1, 0)

    assert grid.map_coord_to_bin_idx(
        (0.75, theta),
    ) is None


# =============================================================================
# Region selection tests
# =============================================================================


def test_get_bins_in_region_non_wrapping_region():
    """
    Test get_bins_in_region without angular wraparound.

    Test that it:
        1) returns all radial bins in region
        2) returns all angular bins in region
        3) returns BinAddress objects
    """
    grid = PolarBinGrid(0.0, 1.0, 4, 8)

    bins = grid.get_bins_in_region(
        (0.1, 0.1),
        (0.6, 0.6),
        False,
    )

    expected = {
        BinAddress(0, 0),
        BinAddress(1, 0),
        BinAddress(2, 0),
    }

    assert bins == expected
    assert all(
        isinstance(bin_address, BinAddress)
        for bin_address in bins
    )


def test_get_bins_in_region_reversed_corners():
    """
    Test get_bins_in_region with reversed corners.

    Test that it:
        1) handles decreasing radial coordinates
        2) handles decreasing angular coordinates
        3) produces identical output
    """
    grid = PolarBinGrid(0.0, 1.0, 4, 8)

    forward_bins = grid.get_bins_in_region(
        (0.1, 0.1),
        (0.6, 0.6),
        False,
    )

    reversed_bins = grid.get_bins_in_region(
        (0.6, 0.6),
        (0.1, 0.1),
        False,
    )

    assert reversed_bins == forward_bins


def test_get_bins_in_region_wraparound_includes_boundary_bins():
    """
    Test get_bins_in_region with angular wraparound.

    Test that it:
        1) includes bins near theta equals zero
        2) includes bins near theta equals 2pi
        3) excludes bins in middle angular interval
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 8)

    bins = grid.get_bins_in_region(
        (0.1, 1.9 * np.pi),
        (0.9, 0.1 * np.pi),
        True,
    )

    assert grid.map_coord_to_bin_idx(
        (0.5, 0.01),
    ) in bins

    assert grid.map_coord_to_bin_idx(
        (0.5, 1.99 * np.pi),
    ) in bins

    assert grid.map_coord_to_bin_idx(
        (0.5, np.pi),
    ) not in bins


def test_get_bins_in_region_wraparound_matching_theta_bins_selects_all_theta():
    """
    Test wrapped selection beginning and ending in same theta bin.

    Test that it:
        1) treats selection as spanning all theta bins
        2) includes all angular bins
        3) includes all requested radial bins
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 8)

    bins = grid.get_bins_in_region(
        (0.1, 0.1),
        (0.9, 0.2),
        True,
    )

    expected = {
        BinAddress(r_idx, theta_idx)
        for r_idx in range(grid.r.n_bins)
        for theta_idx in range(grid.theta.n_bins)
    }

    assert bins == expected


def test_get_bins_in_region_rejects_out_of_bounds_corners():
    """
    Test get_bins_in_region with out-of-bounds corners.

    Test that it:
        1) rejects first corner outside grid
        2) rejects second corner outside grid
        3) raises clear ValueError
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 4)

    with pytest.raises(
        ValueError,
        match="outside of the grid domain",
    ):
        grid.get_bins_in_region(
            (2.0, 0.1),
            (0.5, 0.2),
            False,
        )

    with pytest.raises(
        ValueError,
        match="outside of the grid domain",
    ):
        grid.get_bins_in_region(
            (0.5, 0.1),
            (2.0, 0.2),
            False,
        )


# =============================================================================
# Edge tests
# =============================================================================


def test_list_all_exposed_edges_empty_input_has_no_edges():
    """
    Test list_all_exposed_edges with empty input.

    Test that it:
        1) accepts empty iterable
        2) returns empty list
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 4)

    assert grid.list_all_exposed_edges([]) == []


def test_list_all_exposed_edges_single_bin_has_four_edges():
    """
    Test list_all_exposed_edges with single bin.

    Test that it:
        1) returns outer edge
        2) returns inner edge
        3) returns left edge
        4) returns right edge
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 4)

    edges = grid.list_all_exposed_edges(
        {BinAddress(0, 0)},
    )

    assert len(edges) == 4


def test_list_all_exposed_edges_two_theta_adjacent_bins_share_internal_edge():
    """
    Test adjacent theta bins.

    Test that it:
        1) omits internal edge
        2) preserves exposed edges
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 4)

    edges = grid.list_all_exposed_edges(
        {
            BinAddress(0, 0),
            BinAddress(0, 1),
        }
    )

    assert len(edges) == 6


def test_list_all_exposed_edges_two_radially_adjacent_bins_share_internal_edge():
    """
    Test adjacent radial bins.

    Test that it:
        1) omits internal edge
        2) preserves exposed edges
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 4)

    edges = grid.list_all_exposed_edges(
        {
            BinAddress(0, 0),
            BinAddress(1, 0),
        }
    )

    assert len(edges) == 6


def test_list_all_exposed_edges_theta_wraparound_bins_share_internal_edge():
    """
    Test theta wraparound adjacency.

    Test that it:
        1) treats first and last theta bins as neighbors
        2) omits internal edge
        3) preserves exposed edges
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 4)

    edges = grid.list_all_exposed_edges(
        {
            BinAddress(0, 0),
            BinAddress(0, 3),
        }
    )

    assert len(edges) == 6


# =============================================================================
# Bin edge tests
# =============================================================================


def test_determine_bin_edge_returns_expected_coordinates():
    """
    Test _determine_bin_edge.

    Test that it:
        1) returns expected outer edge
        2) returns expected inner edge
        3) returns expected left edge
        4) returns expected right edge
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 4)

    bin_address = BinAddress(0, 0)

    outer = grid._determine_bin_edge(
        bin_address,
        "outer",
    )

    np.testing.assert_allclose(
        outer,
        ((0.5, 0.0), (0.5, 0.5 * np.pi)),
    )

    inner = grid._determine_bin_edge(
        bin_address,
        "inner",
    )

    np.testing.assert_allclose(
        inner,
        ((0.0, 0.0), (0.0, 0.5 * np.pi)),
    )

    left = grid._determine_bin_edge(
        bin_address,
        "left",
    )

    np.testing.assert_allclose(
        left,
        ((0.0, 0.0), (0.5, 0.0)),
    )

    right = grid._determine_bin_edge(
        bin_address,
        "right",
    )

    np.testing.assert_allclose(
        right,
        ((0.0, 0.5 * np.pi), (0.5, 0.5 * np.pi)),
    )


def test_determine_bin_edge_rejects_unknown_side():
    """
    Test _determine_bin_edge input validation.

    Test that it:
        1) rejects invalid edge type
        2) raises clear ValueError
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 4)

    with pytest.raises(
        ValueError,
        match="Unknown edge type",
    ):
        grid._determine_bin_edge(
            BinAddress(0, 0),
            "diagonal",
        )


# =============================================================================
# Bin area tests
# =============================================================================


def test_calc_bin_area_returns_expected_area_for_inner_radial_bin():
    """
    Test calc_bin_area for inner radial bin.

    Test that it:
        1) uses radial midpoint
        2) uses radial width
        3) uses angular width
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 4)

    area = grid.calc_bin_area(
        BinAddress(0, 0),
    )

    expected = 0.25 * 0.5 * (0.5 * np.pi)

    np.testing.assert_allclose(
        area,
        expected,
    )


def test_calc_bin_area_returns_expected_area_for_outer_radial_bin():
    """
    Test calc_bin_area for outer radial bin.

    Test that it:
        1) uses correct radial shell
        2) uses correct midpoint
        3) returns expected area
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 4)

    area = grid.calc_bin_area(
        BinAddress(1, 0),
    )

    expected = 0.75 * 0.5 * (0.5 * np.pi)

    np.testing.assert_allclose(
        area,
        expected,
    )


def test_calc_bin_area_is_independent_of_theta_index():
    """
    Test calc_bin_area across theta bins.

    Test that it:
        1) gives identical area for same radial shell
        2) depends only on radial index
    """
    grid = PolarBinGrid(0.0, 1.0, 2, 4)

    area0 = grid.calc_bin_area(
        BinAddress(1, 0),
    )

    area3 = grid.calc_bin_area(
        BinAddress(1, 3),
    )

    np.testing.assert_allclose(
        area0,
        area3,
    )


def test_calc_bin_area_handles_nonzero_r_min():
    """
    Test calc_bin_area with nonzero radial minimum.

    Test that it:
        1) includes r_min
        2) uses correct radial midpoint
    """
    grid = PolarBinGrid(1.0, 3.0, 2, 4)

    area = grid.calc_bin_area(
        BinAddress(0, 0),
    )

    expected = 1.5 * 1.0 * (0.5 * np.pi)

    np.testing.assert_allclose(
        area,
        expected,
    )