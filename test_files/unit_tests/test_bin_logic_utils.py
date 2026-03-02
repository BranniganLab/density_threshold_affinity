#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 14:53:51 2026

@author: js2746
"""

import pytest
import numpy as np

from dta.bin_logic.utils import BinAddress, Coordinate, BinEdge, unwrap_theta


# -----------------------
# NamedTuple structural tests
# -----------------------

def test_binaddress_fields_and_tuple_behavior():
    b = BinAddress(r_index=3, theta_index=7)

    # field names / order are part of the public API
    assert BinAddress._fields == ("r_index", "theta_index")

    # tuple-like behavior
    assert tuple(b) == (3, 7)
    assert b[0] == 3
    assert b[1] == 7

    # attribute access
    assert b.r_index == 3
    assert b.theta_index == 7

    # immutability
    with pytest.raises(AttributeError):
        b.r_index = 99  # type: ignore[misc]

    # hashability (useful as dict keys / set elements)
    d = {b: "ok"}
    assert d[BinAddress(3, 7)] == "ok"


def test_coordinate_fields_and_types_are_reasonable():
    c = Coordinate(r_coord=1.5, theta_coord=0.25)

    assert Coordinate._fields == ("r_coord", "theta_coord")
    assert tuple(c) == (1.5, 0.25)
    assert c.r_coord == 1.5
    assert c.theta_coord == 0.25

    with pytest.raises(AttributeError):
        c.theta_coord = 0.0  # type: ignore[misc]

    # hashability sanity
    s = {c}
    assert Coordinate(1.5, 0.25) in s


def test_binedge_fields_and_nesting():
    p1 = Coordinate(1.0, 0.0)
    p2 = Coordinate(2.0, np.pi / 2)
    e = BinEdge(endpoint1=p1, endpoint2=p2)

    assert BinEdge._fields == ("endpoint1", "endpoint2")
    assert e.endpoint1 == p1
    assert e.endpoint2 == p2
    assert e.endpoint1.r_coord == 1.0
    assert e.endpoint2.theta_coord == pytest.approx(np.pi / 2)

    with pytest.raises(AttributeError):
        e.endpoint1 = p2  # type: ignore[misc]


# -----------------------
# unwrap_theta tests
# -----------------------

def test_unwrap_theta_when_previous_is_none_returns_current_unchanged():
    assert unwrap_theta(None, 0.1) == 0.1
    assert unwrap_theta(None, 2 * np.pi - 1e-6) == pytest.approx(2 * np.pi - 1e-6)


def test_unwrap_theta_no_wrap_when_delta_within_pi():
    prev = 1.0
    cur = 1.0 + (np.pi - 1e-9)  # slightly less than +pi
    assert unwrap_theta(prev, cur) == pytest.approx(cur)

    cur2 = 1.0 - (np.pi - 1e-9)  # slightly less than -pi
    assert unwrap_theta(prev, cur2) == pytest.approx(cur2)


def test_unwrap_theta_wraps_negative_when_delta_greater_than_pi():
    # If current_theta jumps forward across 2π boundary relative to previous_theta,
    # delta can appear > pi and we subtract 2π.
    prev = 2 * np.pi - 0.1
    cur = 0.05  # small angle after wrap
    # delta = 0.05 - (2π - 0.1) = -2π + 0.15, which is < -pi, so should add 2π
    # (this case is actually handled by the delta < -pi branch)
    assert unwrap_theta(prev, cur) == pytest.approx(cur + 2 * np.pi)

    # Construct a delta > +pi case directly
    prev2 = 0.1
    cur2 = 0.1 + (np.pi + 1e-6)
    assert unwrap_theta(prev2, cur2) == pytest.approx(cur2 - 2 * np.pi)


def test_unwrap_theta_wraps_positive_when_delta_less_than_minus_pi():
    prev = 0.1 + (np.pi + 1e-6)
    cur = 0.1
    assert unwrap_theta(prev, cur) == pytest.approx(cur + 2 * np.pi)


def test_unwrap_theta_exactly_pi_or_minus_pi_no_adjustment():
    prev = 1.234

    cur_plus_pi = prev + np.pi
    assert unwrap_theta(prev, cur_plus_pi) == pytest.approx(cur_plus_pi)

    cur_minus_pi = prev - np.pi
    assert unwrap_theta(prev, cur_minus_pi) == pytest.approx(cur_minus_pi)


def test_unwrap_theta_idempotent_for_continuous_sequence():
    # If we feed a sequence already continuous, unwrap should not introduce drift.
    thetas = [0.1, 0.2, 0.3, 0.4]
    prev = None
    out = []
    for t in thetas:
        u = unwrap_theta(prev, t)
        out.append(u)
        prev = u
    assert out == pytest.approx(thetas)
