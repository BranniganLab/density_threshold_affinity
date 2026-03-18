#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for lightweight bin-logic value types and angular unwrapping."""

import numpy as np
import pytest

from dta.bin_logic.utils import BinAddress, BinEdge, Coordinate, unwrap_theta


def test_binaddress_exposes_expected_public_tuple_api():
    """Verify that BinAddress preserves its public field order and tuple semantics."""
    b = BinAddress(r_index=3, theta_index=7)

    # Verify the public field names and ordering expected by callers.
    assert BinAddress._fields == ("r_index", "theta_index")

    # Verify tuple-style and attribute-style access expose the same values.
    assert tuple(b) == (3, 7)
    assert b[0] == 3
    assert b[1] == 7
    assert b.r_index == 3
    assert b.theta_index == 7


def test_coordinate_is_immutable_and_hashable():
    """Verify that Coordinate behaves like an immutable, hashable value object."""
    c = Coordinate(r_coord=1.5, theta_coord=0.25)

    # Verify attribute access exposes the stored coordinate values.
    assert c.r_coord == 1.5
    assert c.theta_coord == 0.25

    # Verify the value object cannot be mutated after construction.
    with pytest.raises(AttributeError):
        c.theta_coord = 0.0  # type: ignore[misc]

    # Verify equivalent values can be used interchangeably in hashed containers.
    s = {c}
    assert Coordinate(1.5, 0.25) in s


def test_binedge_stores_coordinate_endpoints():
    """Verify that BinEdge retains the provided Coordinate endpoints unchanged."""
    p1 = Coordinate(1.0, 0.0)
    p2 = Coordinate(2.0, np.pi / 2)
    e = BinEdge(endpoint1=p1, endpoint2=p2)

    # Verify the edge preserves the supplied endpoint objects.
    assert e.endpoint1 == p1
    assert e.endpoint2 == p2

    # Verify nested coordinate data remains accessible through the edge.
    assert e.endpoint1.r_coord == 1.0
    assert e.endpoint2.theta_coord == pytest.approx(np.pi / 2)


def test_unwrap_theta_returns_current_when_no_previous_angle_exists():
    """Verify that unwrap_theta leaves the first sample unchanged."""
    assert unwrap_theta(None, 0.1) == 0.1
    assert unwrap_theta(None, 2 * np.pi - 1e-6) == pytest.approx(2 * np.pi - 1e-6)


def test_unwrap_theta_leaves_small_deltas_unchanged():
    """Verify that unwrap_theta does not adjust values already within pi of the previous angle."""
    prev = 1.0

    # Verify a delta just below +pi is treated as already continuous.
    cur_forward = 1.0 + (np.pi - 1e-9)
    assert unwrap_theta(prev, cur_forward) == pytest.approx(cur_forward)

    # Verify a delta just above -pi is also treated as already continuous.
    cur_backward = 1.0 - (np.pi - 1e-9)
    assert unwrap_theta(prev, cur_backward) == pytest.approx(cur_backward)


def test_unwrap_theta_adjusts_across_wrap_boundary_in_both_directions():
    """Verify that unwrap_theta removes discontinuities larger than pi in either direction."""
    # Verify an angle just after 0 is lifted into continuity with a prior angle near 2π.
    prev_near_2pi = 2 * np.pi - 0.1
    cur_after_wrap = 0.05
    assert unwrap_theta(prev_near_2pi, cur_after_wrap) == pytest.approx(
        cur_after_wrap + 2 * np.pi
    )

    # Verify a forward jump greater than +pi is pulled back by one full revolution.
    prev_small = 0.1
    cur_large_jump = 0.1 + (np.pi + 1e-6)
    assert unwrap_theta(prev_small, cur_large_jump) == pytest.approx(
        cur_large_jump - 2 * np.pi
    )

    # Verify a backward jump less than -pi is pushed forward by one full revolution.
    prev_large = 0.1 + (np.pi + 1e-6)
    cur_small = 0.1
    assert unwrap_theta(prev_large, cur_small) == pytest.approx(cur_small + 2 * np.pi)


def test_unwrap_theta_does_not_adjust_exactly_pi_boundaries():
    """Verify that unwrap_theta treats exact +/-pi deltas as non-wrapping cases."""
    prev = 1.234

    # Verify the +pi boundary does not trigger a wrap correction.
    cur_plus_pi = prev + np.pi
    assert unwrap_theta(prev, cur_plus_pi) == pytest.approx(cur_plus_pi)

    # Verify the -pi boundary does not trigger a wrap correction.
    cur_minus_pi = prev - np.pi
    assert unwrap_theta(prev, cur_minus_pi) == pytest.approx(cur_minus_pi)


def test_unwrap_theta_preserves_an_already_continuous_sequence():
    """Verify that unwrap_theta does not introduce drift into a continuous series."""
    thetas = [0.1, 0.2, 0.3, 0.4]
    prev = None
    out = []

    # Run the sequence incrementally the same way production code would.
    for theta in thetas:
        unwrapped = unwrap_theta(prev, theta)
        out.append(unwrapped)
        prev = unwrapped

    # Verify the result matches the original continuous input.
    assert out == pytest.approx(thetas)
