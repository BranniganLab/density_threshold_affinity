#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 13:55:09 2024.

@author: js2746
"""
import numpy as np
from scipy import constants
import math


def calculate_dG(counts_histogram, n_peak, temperature):
    """
    Calculate the delta G for bulk or site.

    Parameters
    ----------
    counts_histogram : np.ndarray
        The site or bulk counts histogram attribute of a Site.
    n_peak : int
        The mode of the bulk counts histogram.
    temperature : float
        The temperature of your system in K.

    Returns
    -------
    delta_G : float
        The binding affinity for the site or bulk, in kcal/mol.

    """
    minus_RT = -1.0 * temperature * constants.R / 4184.  # 4184 converts J to kcal
    P_unnoc = calculate_P_unnoc(counts_histogram, n_peak)
    delta_G = minus_RT * np.log((1 - P_unnoc) / P_unnoc)
    return delta_G


def calculate_P_unnoc(counts_histogram, n_peak):
    """
    Calculate the probability that the site is unoccupied by ligand. This \
    is defined by the portion of the histogram <= n_peak divided by the \
    total histogram.

    Parameters
    ----------
    counts_histogram : np.ndarray
        The site or bulk counts histogram attribute of a Site.
    n_peak : int
        The mode of the bulk counts histogram.

    Returns
    -------
    P_unnoc : float
        The probability that the site is unoccupied by ligand.

    """
    total_N = np.sum(counts_histogram)
    P_unnoc = np.sum(counts_histogram[:n_peak + 1]) / total_N
    P_occ = np.sum(counts_histogram[n_peak + 1:]) / total_N
    assert math.isclose(P_occ + P_unnoc, 1, abs_tol=0.01), f"Probabilities do not sum to one. Current sum: {P_unnoc + P_occ}"
    return P_unnoc


def calculate_hist_mode(histogram, nonzero=False):
    """
    Calculate the mode of the counts histogram.

    Parameters
    ----------
    histogram : ndarray
        The histogram whose mode you wish to calculate.
    nonzero : boolean
        If True, in the case of mode=0, use the second highest peak instead.\
        This may be necessary when estimating the predicted accessible area.

    Returns
    -------
    mode : int
        The mode of the distribution.

    """
    mode = np.argmax(histogram)
    if len(np.shape(mode)) != 0:
        print(f"Warning: More than one peak identified ({mode}), using first peak ({mode[0]})")
        mode = mode[0]
    if mode == 0:
        if nonzero:
            print("Warning: found an experimental mode of 0 for site, using second highest peak")
            mode = np.argmax(histogram[1:]) + 1
    return mode


def calculate_hist_mean(counts_data):
    """
    Calculate the mean of the counts histogram.

    Parameters
    ----------
    counts_data : ndarray
        The histogram whose mean you wish to calculate.

    Returns
    -------
    mean : float
        The mean of the counts histogram.

    """
    total_N = np.sum(counts_data)
    sum_i = 0
    for i in range(len(counts_data)):
        sum_i += i * counts_data[i]
    mean = sum_i / total_N
    return mean
