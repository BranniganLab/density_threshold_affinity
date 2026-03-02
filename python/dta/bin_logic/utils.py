#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Core geometric value types and angle utilities used throughout DTA.

This module defines lightweight, immutable data structures for representing
polar-grid concepts (bin addresses, coordinates, and bin edges), along with
numerical helpers for working with angular values.

The NamedTuple types in this module serve as stable, hashable value objects
that are safe to use as dictionary keys, set elements, and components of
higher-level domain models.

Functions
---------
unwrap_theta(previous_theta, current_theta)
    Adjust an angular value to maintain continuity across the 2π boundary.
    This is primarily used during drag operations or other incremental
    updates where angular wraparound would otherwise introduce discontinuities.
"""
from typing import NamedTuple
from typing import Union
import numpy as np


class BinAddress(NamedTuple):
    """The r and theta indices that correspond to one bin in a lattice."""

    r_index: int
    theta_index: int


class Coordinate(NamedTuple):
    """The r and theta coordinate pair."""

    r_coord: float
    theta_coord: float


class BinEdge(NamedTuple):
    """
    Geometric description of a single visible bin edge in polar coordinates.

    Attributes
    ----------
    endpoint1, endpoint2 : Coordinate
        The coordinates of both edge endpoints.
    """

    endpoint1: Coordinate
    endpoint2: Coordinate


def unwrap_theta(previous_theta, current_theta):
    """
    Unwrap a circular angle measurement to produce a continuous angular sequence.

    This function resolves the discontinuity at 0 / 2π by adjusting the current
    angle by ±2π when necessary so that the change relative to the previous angle
    is minimal.

    Parameters
    ----------
    previous_theta : float or None
        The previous angle in the sequence, expressed in radians. If None,
        no unwrapping is performed and current_theta is returned unchanged.
    current_theta : float
        The current angle measurement in radians, typically constrained to
        the interval [0, 2π).

    Returns
    -------
    float
        The unwrapped angle in radians, continuous with previous_theta.
    """
    if previous_theta is None:
        return current_theta

    delta = current_theta - previous_theta

    if delta > np.pi:
        return current_theta - 2 * np.pi
    if delta < -np.pi:
        return current_theta + 2 * np.pi

    return current_theta


CoordinateLike = Union[Coordinate, tuple[float, float]]
BinAddressLike = Union[BinAddress, tuple[int, int]]

def as_coordinate(x: CoordinateLike) -> Coordinate:
    """Coerce (r, theta) tuples into Coordinate; pass through Coordinate."""
    if isinstance(x, Coordinate):
        return x
    r, theta = x
    return Coordinate(r_coord=float(r), theta_coord=float(theta))

def as_bin_address(bin_: BinAddressLike) -> BinAddress:
    """Normalize user input to a BinAddress."""
    if isinstance(bin_, BinAddress):
        return bin_
    r, t = bin_
    return BinAddress(int(r), int(t))
