#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Angular utility classes and functions for polar coordinate systems.

This module provides small, reusable helpers for working with angular
quantities (e.g., wrapping, unwrapping, and normalization) in radians.

Functions here are framework-agnostic and shared by core geometry,
GUI interaction logic, and analysis code.
"""
from typing import NamedTuple
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
