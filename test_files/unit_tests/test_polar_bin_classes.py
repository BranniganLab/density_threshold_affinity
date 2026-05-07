#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 16:01:23 2026

@author: js2746
"""

import numpy as np
from dta.bin_logic import PolarBinGrid


def test_map_coord_to_bin_idx_basic_and_wrap():
    grid = PolarBinGrid(0, 1, 2, 4)

    r = 0.25
    theta = 0.25 * np.pi

    assert grid.map_coord_to_bin_idx((r, theta)) == (0, 0)
    assert grid.map_coord_to_bin_idx((r, theta + 2.0 * np.pi)) == (0, 0)

    # outside radius
    assert grid.map_coord_to_bin_idx((2.0, theta)) is None

    # exactly on outer edge is outside
    assert grid.map_coord_to_bin_idx((1.0, theta)) is None


def test_map_coord_to_bin_idx_boundaries_theta_and_r():
    grid = PolarBinGrid(0, 2, 2, 2)

    # theta exactly on interior edge goes to the bin on the "right"
    assert grid.map_coord_to_bin_idx((0.5, np.pi)) == (0, 1)

    # theta=0.0 maps to first bin
    assert grid.map_coord_to_bin_idx((0.5, 0.0)) == (0, 0)

    # theta at the last *grid edge* is outside (since this grid does not span 2π)
    assert grid.map_coord_to_bin_idx((0.5, 2.0 * np.pi)) is None

    # actual angular wraparound happens at 2π
    assert grid.map_coord_to_bin_idx((0.5, 2.0 * np.pi + 1)) == (0, 0)

    # r exactly on interior edge goes to bin on the "right"
    assert grid.map_coord_to_bin_idx((1.0, 0.5)) == (1, 0)

    # r at outer edge is out of range
    assert grid.map_coord_to_bin_idx((2.0, 0.5)) is None


def test_bins_in_region_nonempty_and_indices_valid():
    grid = PolarBinGrid(0, 1, 4, 8)

    bins = grid.bins_in_region((0.1, 0.1), (0.6, 0.6))

    assert len(bins) > 0
    assert all(0 <= ri < grid.n_r and 0 <= ti < grid.n_theta for (ri, ti) in bins)


def test_bins_in_region_wraparound_includes_zero_angle_bin():
    grid = PolarBinGrid(0, 1, 2, 8)

    # Cross 2pi boundary: near 2pi down to small angle
    bins = grid.bins_in_region((0.0, 1.9 * np.pi), (1.0, 0.1 * np.pi))
    assert len(bins) > 0

    # Must include a bin near theta ~ 0
    assert grid.map_coord_to_bin_idx((0.5, 0.01)) in bins


def test_exposed_edges_single_bin_has_four():
    grid = PolarBinGrid(0, 1, 2, 4)
    edges = grid.exposed_edges({(0, 0)})
    assert len(edges) == 4


def test_exposed_edges_two_adjacent_bins_share_internal_edge():
    grid = PolarBinGrid(0, 1, 2, 4)

    # adjacent in theta at same r
    edges = grid.exposed_edges({(0, 0), (0, 1)})

    # Two squares share one internal edge => 4 + 4 - 2 = 6 exposed edges
    assert len(edges) == 6
