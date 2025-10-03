#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 13:55:09 2024.

@author: js2746
"""
import numpy as np
from scipy import constants
import math
from pathlib import Path
from DTA.density import parse_tcl_dat_file


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


def load_inclusion_helices(path):
    """
    Get helix locations from the given path.

    Parameters
    ----------
    path  :  str or Path
        The path to the directory containing the helix coordinate files.

    Returns
    -------
    helix_list  :  list
        List of upper and lower helix coordinates.
    """
    if isinstance(path, str):
        path = Path(path)
    try:
        helices_lwr = np.loadtxt(path.joinpath("Protein_coords_lwr.dat"))
        helices_upr = np.loadtxt(path.joinpath("Protein_coords_upr.dat"))
    except FileNotFoundError:
        helices_lwr = None
        helices_upr = None
        print("Protein coordinates not found")

    return [helices_upr, helices_lwr]


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


def load_replica_counts(root_path, replicas_list, system_name, leaflet_id):
    """
    Load the counts from all replicas of a single system. Return them as a list.

    Parameters
    ----------
    root_path : str or Path
        The path to directory containing your replica subdirectories.
    replicas_list : list
        The list of replica subdirectory names.
    system_name : str
        The file stem for your PolarDensityBin outputs.
    leaflet_id : int
        Which leaflet your Site is in. 1=outer and 2=inner.

    Returns
    -------
    list of ndarrays.

    """
    assert leaflet_id in [1, 2], "leaflet_id must be 1 (outer) or 2 (inner)."
    assert isinstance(root_path, (str, Path)), "root_path must be a str or Path."
    if isinstance(root_path, str):
        root_path = Path(root_path)
    assert root_path.exists(), f"could not find root_path {root_path}"
    assert isinstance(replicas_list, list), "replicas_list must be a list."
    assert (len(replicas_list) > 1), "Less than 2 replicas found."
    replica_counts_list = []
    leaflet = {1: "upp", 2: "low"}
    for rep in replicas_list:
        fname = root_path.joinpath(rep, f"{system_name}.{leaflet[leaflet_id]}.dat")
        assert fname.is_file(), f"could not find file {fname}"
        counts, grid_dims, system_info = parse_tcl_dat_file(fname, bulk=False)
        replica_counts_list.append(counts)
    return replica_counts_list
