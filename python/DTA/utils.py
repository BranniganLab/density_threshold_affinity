#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 13:55:09 2024.

@author: js2746
"""
import math
import warnings
from pathlib import Path
import numpy as np
from scipy import constants


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
    P_unocc = calculate_P_unocc(counts_histogram, n_peak)
    delta_G = minus_RT * np.log((1 - P_unocc) / P_unocc)
    return delta_G


def calculate_P_unocc(counts_histogram, n_peak):
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
    P_unocc : float
        The probability that the site is unoccupied by ligand.

    """
    total_N = np.sum(counts_histogram)
    P_unocc = np.sum(counts_histogram[:n_peak + 1]) / total_N
    P_occ = np.sum(counts_histogram[n_peak + 1:]) / total_N
    assert math.isclose(P_occ + P_unocc, 1, abs_tol=0.01), f"Probabilities do not sum to one. Current sum: {P_unocc + P_occ}"
    return P_unocc


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
    for index, value in enumerate(counts_data):
        sum_i += index * value
    mean = sum_i / total_N
    return mean


def load_inclusion_coordinates(directory):
    """
    Get backbone coordinates from the given directory.

    Parameters
    ----------
    directory  :  str or Path
        The path to the directory containing the helix coordinate files.

    Returns
    -------
    backbone_coms  :  2-list of 1D numpy ndarrays
        List of outer and inner leaflet backbone coordinates, alternating r and theta.
    """
    path = validate_path(directory)
    backbone_com_upr = None
    backbone_com_lwr = None
    fails = 0
    for leaflet in ["upr", "lwr"]:
        fname = path.joinpath(f"Protein_coords_{leaflet}.dat")
        try:
            with open(fname, 'r', encoding="utf-8") as file:
                coords = file.readlines()
        except FileNotFoundError as err:
            # Sometimes user will only want coordinates from one leaflet. Only
            # error if there are no coordinates from both leaflets.
            coords = None
            warnings.warn(f"No coordinates found in {leaflet} leaflet.")
            fails += 1
            if fails == 2:
                raise FileNotFoundError(f"Could not find protein coordinate files in {path}") from err
        coords_list = remove_vals_that_match_substring(coords, "/")
        if leaflet == "upr":
            backbone_com_upr = coords_list
        else:
            backbone_com_lwr = coords_list
    return [backbone_com_upr, backbone_com_lwr]


def remove_vals_that_match_substring(input_text, substring):
    """
    Iterate through lines and remove anything that contains substring.

    Parameters
    ----------
    input_text : list
        The text that was read-in from file with readlines().
    substring : str
        Filter out any values that contain substring.

    Returns
    -------
    list of lists
        The filtered text.

    """
    if input_text is None:
        return []
    output_text = []
    for row in input_text:
        row = row.strip()
        if row:
            if row[0] == "#":
                continue
            temp = []
            for value in row.split(' '):
                if substring in value:
                    continue
                temp.append(value)
            output_text.append(temp)
    return output_text


def aggregate_site_counts_histograms(site_list):
    """
    Cycle through all the sites and add their counts histograms together.

    Parameters
    ----------
    site_list : list
        List of Sites.

    Returns
    -------
    counts : ndarray
        1D numpy array of histogrammed bead counts.

    """
    hist_lengths = []
    for site in site_list:
        assert site.site_counts_histogram is not None, "One or more sites do not have counts associated. Please use update_counts_histogram() and try again."
        hist_length = site.site_counts_histogram.shape[0]
        hist_lengths.append(hist_length)
    max_len = max(hist_lengths)
    counts = np.zeros(max_len)
    for site in site_list:
        counts_to_add = site.site_counts_histogram.copy()
        if counts_to_add.shape[0] < max_len:
            padding = max_len - counts_to_add.shape[0]
            counts_to_add = np.pad(counts_to_add, (0, padding), mode='constant', constant_values=0)
        counts += counts_to_add
    return counts


def check_bulk_counts_histogram(site_list):
    """
    Cycle through each Site and make sure the bulk_counts_histograms all match.\
    Return one of them.

    Parameters
    ----------
    site_list : list
        List of Sites.

    Returns
    -------
    bulk : ndarray
        1D numpy array of histogrammed bead counts from the bulk distribution.

    """
    first_site = site_list[0]
    bulk = first_site.bulk_counts_histogram.copy()
    for site in site_list[1:]:
        assert bulk.all() == site.bulk_counts_histogram.all(), "One or more sites have different bulk histograms. This shouldn't be possible."
    return bulk


def validate_path(path, file=False):
    """
    Make sure that a path can be found. If file=True, make sure that it is \
        recognized as a file. If path is a string, convert it to a Pathlib Path.

    Parameters
    ----------
    path : str or Path
        The path you would like to validate.
    file : bool, optional
        If True, check to see if the file is recognized. The default is False.

    Raises
    ------
    FileNotFoundError
    TypeError

    Returns
    -------
    path : Path
        Returns the path as a Path object.

    """
    if isinstance(path, str):
        path = Path(path)
    path = path.resolve()
    if not isinstance(path, Path):
        raise TypeError(f"path must be a string or a Pathlib Path. Object of type {type(path)} was supplied.")
    if not path.exists():
        raise FileNotFoundError(f"The specified path '{path}' could not be found. Please ensure the path is correct.")
    if file:
        if not path.is_file():
            raise FileNotFoundError(f"The specified file '{path}' is not recognized as a file. Please ensure the path is correct.")
    return path


def theta_in_bin(t0, t1, t_low, t_high):
    dt = (t1 - t0) % (2 * np.pi)
    theta_val = t0
    if dt > np.pi:
        dt = (t0 - t1) % (2 * np.pi)
        theta_val = t1
    condition1 = (t_low - theta_val) % (2 * np.pi)
    condition2 = (t_high - theta_val) % (2 * np.pi)
    condition3 = (theta_val - t_low) % (2 * np.pi)
    return ((condition1 < dt) or (condition2 < dt) or (condition3 < (t_high - t_low)))
