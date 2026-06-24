#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 09:44:53 2024.

@author: js2746
"""
from pathlib import Path
from collections import namedtuple
import numpy as np
from dta.utils import validate_path
from dta.bin_logic import PolarBinGrid


SysInfo = namedtuple('SysInfo', ['NL', 'NB', 'BoxArea', 'ExpBeadDensity', 'DrDtheta'])
Dimensions = namedtuple('Dimensions', ['dr', 'Nr', 'dtheta', 'Ntheta', 'Nframes'])


def parse_tcl_dat_file(filepath, bulk):
    """
    Extract bin counts and bin dimensions from a TCL output .dat file.

    Parameters
    ----------
    filepath : Path
        The path to the TCL .dat or .out output file.
    bulk : boolean
        If True, assume output is from do_get_counts. If False, assume \
        output is from polarDensityBin.

    Returns
    -------
    counts : ndarray
        If bulk is True, return a 1D ndarray. If bulk is False, return a 3D \
        ndarray of integer counts. Dimension 0 is time, dimension 1 is r and \
        dimension 2 is theta. IE the count for the 4th radial bin and 12th \
        theta bin on the 33rd frame of the trajectory would be counts[32, 3, 11].
    dimensions : namedtuple or None
        If bulk is True, return None. If bulk is False, return the bin \
        dimensions in r and theta as well as the number of frames inside a \
        namedtuple.
    system_info : namedtuple or None
        The information taken from the header row of the .dat file specified in\
        filepath. If bulk is True, returns None.

    """
    filepath = validate_path(filepath, file=True)
    assert (filepath.suffixes[-1] == '.dat') or (filepath.suffixes[-1] == '.out'), "You must provide the .dat or .out file output from VMD."
    if bulk:
        return np.loadtxt(filepath).astype(int).flatten(), None, None
    unrolled_data = np.loadtxt(filepath, skiprows=1)
    system_info = _parse_system_info(np.loadtxt(filepath, comments=None, max_rows=1, delimiter=',', dtype=str))
    grid = _create_grid_from_dat_file(unrolled_data)
    counts = _package_counts(unrolled_data, grid).squeeze()
    return counts, grid, system_info


def valid_Dimensions(list_of_Dimensions_objs):
    """
    Make sure that all Dimensions attributes do not vary across objects \
        (not counting Nframes) in a list of Dimensions objects. If attributes \
        do not vary across objects, return True. Else, return False.

    Parameters
    ----------
    list_of_Dimensions_objs : list
        List containing multiple Dimensions objects.

    Returns
    -------
    Boolean

    """
    if not isinstance(list_of_Dimensions_objs, list):
        raise TypeError(f"grid_dims must be a list, not a {type(list_of_Dimensions_objs)}.")
    if not all(isinstance(item, Dimensions) for item in list_of_Dimensions_objs):
        raise TypeError("grid_dims_list must contain Dimensions namedtuples.")
    if len(list_of_Dimensions_objs) == 1:
        return True
    dr, Nr, dtheta, Ntheta, _ = list_of_Dimensions_objs[0]
    for obj in list_of_Dimensions_objs:
        for index, val in enumerate([dr, Nr, dtheta, Ntheta]):
            if obj[index] != val:
                return False
    return True


def load_replica_counts(root_path, replicas_list, system_name, leaflet_id, avg=False):
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
    avg : bool
        If True, load the average density file. If False (default), load the \
        individual frame counts.

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
        if avg:
            fname = root_path.joinpath(rep, f"{system_name}.{leaflet[leaflet_id]}.avg.dat")
        else:
            fname = root_path.joinpath(rep, f"{system_name}.{leaflet[leaflet_id]}.dat")
        assert fname.is_file(), f"could not find file {fname}"
        counts, _, _ = parse_tcl_dat_file(fname, bulk=False)
        replica_counts_list.append(counts)
    return replica_counts_list


def calculate_density(avg_counts, grid):
    """
    Calculate the average bead density in each bin.

    Parameters
    ----------
    avg_counts : ndarray
        2D array of average bead counts per bin in simulation. [r,theta] format.
    grid : PolarBinGrid
        The PolarBinGrid object corresponding to this system.

    Returns
    -------
    density : ndarray
        The average bead density in each bin.

    """
    assert isinstance(avg_counts, np.ndarray), "avg_counts must be an ndarray"
    assert len(avg_counts.shape) == 2, "avg_counts must be a 2D array."
    density = avg_counts / grid.bin_areas
    return density


def calculate_density_enrichment(density, expected_density):
    """
    Divide the average bead densities in each bin by the expected bead density \
    to calculate an enrichment score.

    Parameters
    ----------
    density : ndarray
        2D array of average bead density values.
    expected_density : float
        The expected bead density; output from polarDensityBin.

    Returns
    -------
    ndarray
        Array of same shape as density containing density enrichment scores.

    """
    assert isinstance(density, np.ndarray), "density must be a numpy array."
    assert len(density.shape) == 2, "density must be a 2d array."
    assert isinstance(expected_density, (float, int)), "expected_density must be a scalar."
    return density / expected_density


def aggregate_density_enrichment_scores(file_paths):
    """
    Collect the average counts, bin sizing, and expected density from each \
        system, calulate average density enrichment for each system, and calculate\
        average density enrichment accross systems. Also, return the Dimensions \
        object that needs to be used for plotting the heatmaps.

    Parameters
    ----------
    file_paths : list
        The list of paths to all the average counts files you want to calculate\
        density enrichment for.

    Returns
    -------
    all_reps_avg : numpy ndarray
        The density enrichment scores, averaged together in a 2D array.
    grid_dims_final : namedtuple
        The Dimensions namedtuple that contains information about lattice sizes.

    """
    replica_enrichments_list = []
    replica_dims_list = []
    for rep_path in file_paths:
        rep_path = validate_path(rep_path, file=True)
        counts, grid_dims, system_info = parse_tcl_dat_file(rep_path, bulk=False)
        replica_dims_list.append(grid_dims)
        density_enrichment = calculate_density_enrichment(calculate_density(counts, grid_dims), system_info.ExpBeadDensity)
        replica_enrichments_list.append(density_enrichment)

    if not valid_Dimensions(replica_dims_list):
        raise ValueError("Not all Dimensions attributes match.")

    all_reps_avg = np.nanmean(np.stack(tuple(replica_enrichments_list), axis=0), axis=0)
    grid_dims_final = replica_dims_list[0]

    return all_reps_avg, grid_dims_final


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  UTIL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%% #


def _parse_system_info(dat_file_header):
    """
    Turn the list of strings from the .dat file header into a namedtuple with \
    useable float values.

    Parameters
    ----------
    dat_file_header : list
        The header line of the .dat file output by polarDensityBin, saved as a \
        list of strings.

    Returns
    -------
    sysInfo : namedtuple
        Namedtuple containing system information output by polarDensityBin.

    """
    # this is done for backwards-compatibility with a discontinued feature
    if len(dat_file_header) == 6:
        NL, NB, BoxArea, ExpBeadDensity, _, DrDtheta = dat_file_header
    elif len(dat_file_header) == 5:
        NL, NB, BoxArea, ExpBeadDensity, DrDtheta = dat_file_header
    else:
        raise ValueError("The header in your .dat file has an unexpected number of entries.")
    NL = _isolate_number_from_header_string(NL)
    NB = _isolate_number_from_header_string(NB)
    BoxArea = _isolate_number_from_header_string(BoxArea)
    ExpBeadDensity = _isolate_number_from_header_string(ExpBeadDensity)
    DrDtheta = _isolate_number_from_header_string(DrDtheta)
    sysInfo = SysInfo(NL, NB, BoxArea, ExpBeadDensity, DrDtheta)
    return sysInfo


def _isolate_number_from_header_string(string):
    """
    Return the number following the colon in the input string. \
    Example: 'sample header info : 185.72 sample units' should return 185.72.

    Parameters
    ----------
    string : str
        A string from the header of the .dat file output by polarDensityBin.

    Returns
    -------
    float
        The number following the colon in the input string.

    """
    right_side = string.split(':')[1].strip()
    if "/" in right_side:
        # Expected bead density has a division symbol in it
        return float(right_side.split('/')[0])
    return float(right_side.split()[0])


def _create_grid_from_dat_file(dat_file_contents):
    """
    Return the length and number of bins in each dimension, as well as the \
    number of frames in the .dat file.

    Parameters
    ----------
    dat_file_contents : ndarray
        The data from the .dat file output by polarDensityBin, minus its header row.

    Returns
    -------
    grid : PolarBinGrid
        Contains r_min, r_max, n_r, and n_theta.

    """
    r_min = dat_file_contents[0, 0]
    d_r = dat_file_contents[1, 0] - r_min
    r_max = np.max(dat_file_contents[1, :])
    n_r = int(round((r_max - r_min) / d_r))
    d_theta_degrees = dat_file_contents[0, 2]
    n_theta = int(round(360 / d_theta_degrees))
    grid = PolarBinGrid(r_min, r_max, n_r, n_theta)
    return grid


def _package_counts(unrolled_data, grid):
    """
    Package the unrolled data into a 3d array in [time, r, theta] format.

    Parameters
    ----------
    unrolled_data : ndarray
        2d array output from polarDensityBin, minus header row.
    grid : PolarBinGrid
        Contains bin dimensions and number of frames.

    Returns
    -------
    counts : ndarray
        3d array of counts in [time, r, theta] format.

    """
    nframes = int(round(unrolled_data.shape[0] / grid.n_r))

    # chop off the first few columns
    unrolled_counts = unrolled_data[:, 3:]

    # 'sideways' because it is in [r, time, theta] format at first
    sideways_counts = np.zeros((grid.n_r, nframes, grid.n_theta))
    for i in range(grid.n_r):
        sideways_counts[i, :, :] = unrolled_counts[(nframes * i):(nframes * (i + 1)), :]

    # swap axes to put it in [time, r, theta] format
    counts = np.swapaxes(sideways_counts, 0, 1)

    return counts
