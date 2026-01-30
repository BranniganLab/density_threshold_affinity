#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 17:34:25 2026

@author: js2746
"""

import numpy as np
import pytest

from DTA.polar_bin_classes import PolarBinGrid


TWO_PI = 2.0 * np.pi


@pytest.fixture()
def grid():
    # Edges aren't used by bin_in_theta_arc, but PolarBinGrid needs valid arrays.
    return PolarBinGrid(theta_edges=np.linspace(0.0, TWO_PI, 5), r_edges=np.linspace(0.0, 1.0, 3))


# -----------------------------------------------------------------------------
# Branch: Full circle selects everything
# -----------------------------------------------------------------------------
@pytest.mark.parametrize("dtheta", [TWO_PI, -TWO_PI, TWO_PI + 1e-6, -(TWO_PI + 1e-6)])
def test_full_circle_selects_everything(grid, dtheta):
    theta_start = 1.234
    theta_end = theta_start + dtheta
    # Any bin interval should be selected
    assert grid.bin_in_theta_arc(theta_start, theta_end, 0.1, 0.2) is True
    assert grid.bin_in_theta_arc(theta_start, theta_end, 5.0, 6.0) is True  # arbitrary values


# -----------------------------------------------------------------------------
# Branch: bin_end normalization (<= 2π uses raw value; > 2π uses modulo)
# We specifically verify that bin_end == 2π is NOT modulo'd to 0.
# -----------------------------------------------------------------------------
def test_bin_end_exactly_two_pi_not_modulo(grid):
    # CCW no-wrap arc: [0, 1)
    theta_start = 0.0
    theta_end = 1.0

    # Bin that spans [~2π-0.1, 2π] should NOT intersect [0,1) if bin_end is treated as 2π.
    # If bin_end were modulo'd, it'd become 0 and this would (incorrectly) appear to intersect.
    bin_start = TWO_PI - 0.1
    bin_end = TWO_PI
    assert grid.bin_in_theta_arc(theta_start, theta_end, bin_start, bin_end) is False


def test_bin_end_greater_than_two_pi_is_moduloed(grid):
    # CCW wrap arc: start at 2π-0.2, end at (start + 0.4) wraps; covers [2π-0.2, 2π) U [0, 0.2)
    theta_start = TWO_PI - 0.2
    theta_end = theta_start + 0.4  # wraps

    # Put a bin_end > 2π. After modulo it should behave like a small bin near 0.
    bin_start = 0.05
    bin_end = TWO_PI + 0.15  # > 2π => modulo -> 0.15
    assert grid.bin_in_theta_arc(theta_start, theta_end, bin_start, bin_end) is True


# -----------------------------------------------------------------------------
# Branch: dtheta >= 0 (CCW) and theta_end <= 2π (no wrap)
# Return: not (bin_end_n <= theta_start or bin_start_n >= theta_end)
# Exercise both True and False outcomes.
# -----------------------------------------------------------------------------
def test_ccw_no_wrap_intersects_true(grid):
    # Arc: [1.0, 2.0)
    theta_start = 1.0
    theta_end = 2.0

    # Bin: [1.5, 1.8) intersects
    assert grid.bin_in_theta_arc(theta_start, theta_end, 1.5, 1.8) is True


def test_ccw_no_wrap_intersects_false_bin_before_arc(grid):
    # Arc: [1.0, 2.0)
    theta_start = 1.0
    theta_end = 2.0

    # Bin ends at/ before arc start -> no intersection
    assert grid.bin_in_theta_arc(theta_start, theta_end, 0.2, 1.0) is False
    assert grid.bin_in_theta_arc(theta_start, theta_end, 0.2, 0.9) is False


def test_ccw_no_wrap_intersects_false_bin_after_arc(grid):
    # Arc: [1.0, 2.0)
    theta_start = 1.0
    theta_end = 2.0

    # Bin starts at/after arc end -> no intersection
    assert grid.bin_in_theta_arc(theta_start, theta_end, 2.0, 2.2) is False
    assert grid.bin_in_theta_arc(theta_start, theta_end, 2.1, 2.3) is False


# -----------------------------------------------------------------------------
# Branch: dtheta >= 0 (CCW) and theta_end > 2π (wrap past 2π)
# Return: not ((bin_end_n <= theta_start) and (bin_start_n >= theta_end - 2π))
# Exercise both True and False outcomes.
# -----------------------------------------------------------------------------
def test_ccw_wrap_intersects_true_near_zero_segment(grid):
    # Arc wraps: start 2π-0.2, end 2π+0.2 => covers [2π-0.2, 2π) U [0, 0.2)
    theta_start = TWO_PI - 0.2
    theta_end = TWO_PI + 0.2

    # Bin in [0.05, 0.10) should intersect
    assert grid.bin_in_theta_arc(theta_start, theta_end, 0.05, 0.10) is True


def test_ccw_wrap_intersects_true_near_two_pi_segment(grid):
    theta_start = TWO_PI - 0.2
    theta_end = TWO_PI + 0.2

    # Bin in [2π-0.15, 2π-0.05) should intersect
    assert grid.bin_in_theta_arc(theta_start, theta_end, TWO_PI - 0.15, TWO_PI - 0.05) is True


def test_ccw_wrap_intersects_false_in_gap_between_segments(grid):
    # Arc wraps: [2π-0.2, 2π) U [0, 0.2)
    theta_start = TWO_PI - 0.2
    theta_end = TWO_PI + 0.2
    gap_start = 0.3
    gap_end = 0.5

    # Bin entirely in the "gap" [0.2, 2π-0.2) should NOT intersect
    assert grid.bin_in_theta_arc(theta_start, theta_end, gap_start, gap_end) is False


# -----------------------------------------------------------------------------
# Branch: dtheta < 0 (CW) and theta_end >= 0 (no wrap)
# Return: not (bin_end_n <= theta_end or bin_start_n >= theta_start)
# Exercise both True and False outcomes.
# Note: CW arc is directed from start down to end; geometrically it covers [theta_end, theta_start]
# -----------------------------------------------------------------------------
def test_cw_no_wrap_intersects_true(grid):
    # CW arc from 2.0 down to 1.0 => covers [1.0, 2.0]
    theta_start = 2.0
    theta_end = 1.0

    # Bin inside interval
    assert grid.bin_in_theta_arc(theta_start, theta_end, 1.5, 1.8) is True


def test_cw_no_wrap_intersects_false_bin_before_interval(grid):
    theta_start = 2.0
    theta_end = 1.0

    # Bin entirely below theta_end
    assert grid.bin_in_theta_arc(theta_start, theta_end, 0.2, 0.9) is False
    assert grid.bin_in_theta_arc(theta_start, theta_end, 0.2, 1.0) is False


def test_cw_no_wrap_intersects_false_bin_after_interval(grid):
    theta_start = 2.0
    theta_end = 1.0

    # Bin starts at/after theta_start
    assert grid.bin_in_theta_arc(theta_start, theta_end, 2.0, 2.2) is False
    assert grid.bin_in_theta_arc(theta_start, theta_end, 2.1, 2.3) is False


# -----------------------------------------------------------------------------
# Branch: dtheta < 0 (CW) and theta_end < 0 (wrap below 0)
# Return: not ((bin_end_n <= theta_end + 2π) and (bin_start_n >= theta_start))
# Exercise both True and False outcomes.
# -----------------------------------------------------------------------------
def test_cw_wrap_intersects_true_near_two_pi_segment(grid):
    # CW wrap: start at 0.2, end at -0.2 => covers [0,0.2] U [2π-0.2, 2π)
    theta_start = 0.2
    theta_end = -0.2

    # Bin near 2π should intersect
    assert grid.bin_in_theta_arc(theta_start, theta_end, TWO_PI - 0.15, TWO_PI - 0.05) is True


def test_cw_wrap_intersects_true_near_zero_segment(grid):
    theta_start = 0.2
    theta_end = -0.2

    # Bin near 0 should intersect
    assert grid.bin_in_theta_arc(theta_start, theta_end, 0.05, 0.10) is True


def test_cw_wrap_intersects_false_in_gap_between_segments(grid):
    # Gap is (0.2, 2π-0.2)
    theta_start = 0.2
    theta_end = -0.2

    # Bin in the gap should not intersect
    assert grid.bin_in_theta_arc(theta_start, theta_end, 1.0, 1.2) is False
