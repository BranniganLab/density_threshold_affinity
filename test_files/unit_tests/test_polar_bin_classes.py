#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 16:01:23 2026

@author: js2746
"""

import numpy as np
import pytest

from DTA.polar_bin_classes import PolarBinGrid, BinSelectionModel, PolarBinRenderer, BinEdge


def test_bin_at_basic_and_wrap():
    grid = PolarBinGrid(
        theta_edges=np.linspace(0.0, 2.0 * np.pi, 5),  # 4 bins
        r_edges=np.linspace(0.0, 1.0, 3),              # 2 bins
    )

    r = 0.25
    theta = 0.25 * np.pi

    assert grid.bin_at(r, theta) == (0, 0)
    assert grid.bin_at(r, theta + 2.0 * np.pi) == (0, 0)

    # outside radius
    assert grid.bin_at(2.0, theta) is None

    # exactly on outer edge is outside
    assert grid.bin_at(1.0, theta) is None


def test_bin_at_boundaries_theta_and_r():
    theta_edges = np.array([0.0, 1.0, 2.0])  # 2 theta bins: [0,1), [1,2)
    r_edges = np.array([0.0, 1.0, 2.0])      # 2 r bins: [0,1), [1,2)
    grid = PolarBinGrid(theta_edges, r_edges)

    # theta exactly on interior edge goes to the bin on the "right"
    assert grid.bin_at(0.5, 1.0) == (0, 1)

    # theta=0.0 maps to first bin
    assert grid.bin_at(0.5, 0.0) == (0, 0)

    # theta at the last *grid edge* is outside (since this grid does not span 2π)
    assert grid.bin_at(0.5, 2.0) is None

    # actual angular wraparound happens at 2π
    assert grid.bin_at(0.5, 2.0 * np.pi) == (0, 0)

    # r exactly on interior edge goes to bin on the "right"
    assert grid.bin_at(1.0, 0.5) == (1, 0)

    # r at outer edge is out of range
    assert grid.bin_at(2.0, 0.5) is None


def test_bins_in_region_nonempty_and_indices_valid():
    theta_edges = np.linspace(0, 2*np.pi, 9)  # 8 bins
    r_edges = np.linspace(0, 1, 5)           # 4 bins
    grid = PolarBinGrid(theta_edges, r_edges)

    bins = set(grid.bins_in_region(0.1, 0.1, 0.6, 0.6))

    assert len(bins) > 0
    assert all(0 <= ri < grid.n_r and 0 <= ti < grid.n_t for (ri, ti) in bins)


def test_bins_in_region_wraparound_includes_zero_angle_bin():
    theta_edges = np.linspace(0, 2*np.pi, 9)  # 8 bins
    r_edges = np.linspace(0, 1, 3)            # 2 bins
    grid = PolarBinGrid(theta_edges, r_edges)

    # Cross 2pi boundary: near 2pi down to small angle
    bins = set(grid.bins_in_region(0.0, 1.9*np.pi, 1.0, 0.1*np.pi))
    assert len(bins) > 0

    # Must include a bin near theta ~ 0
    assert grid.bin_at(0.5, 0.01) in bins


def test_exposed_edges_single_bin_has_four():
    grid = PolarBinGrid(
        theta_edges=np.linspace(0, 2*np.pi, 5),
        r_edges=np.linspace(0, 1, 3),
    )
    edges = grid.exposed_edges({(0, 0)})
    assert len(edges) == 4


def test_exposed_edges_two_adjacent_bins_share_internal_edge():
    grid = PolarBinGrid(
        theta_edges=np.linspace(0, 2*np.pi, 5),  # 4 theta bins
        r_edges=np.linspace(0, 1, 3),            # 2 r bins
    )

    # adjacent in theta at same r
    edges = grid.exposed_edges({(0, 0), (0, 1)})

    # Two squares share one internal edge => 4 + 4 - 2 = 6 exposed edges
    assert len(edges) == 6


def test_renderer_draw_edges_creates_artists():
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    renderer = PolarBinRenderer(ax, plot_kwargs={"color": "k", "lw": 1})

    edges = [
        BinEdge((0.0, 1.0), (0.0, 0.0)),
        BinEdge((1.0, 1.0), (0.0, 0.5)),
    ]

    artists = renderer.draw_edges(edges)
    assert len(artists) == 2
    assert all(a.axes is ax for a in artists)


def test_selection_model_ops():
    m = BinSelectionModel()

    m.set({(0, 0), (1, 2)})
    assert m.snapshot() == frozenset({(0, 0), (1, 2)})

    m.add({(2, 3)})
    assert (2, 3) in m.bins()

    m.remove({(0, 0)})
    assert (0, 0) not in m.bins()

    m.clear()
    assert m.bins() == set()
