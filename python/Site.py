#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 13:44:55 2024.

@author: js2746
"""

import numpy as np
import warnings
import math
from pathlib import Path
from scipy import constants
from collections import namedtuple
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
import matplotlib.gridspec as gridspec

Dimensions = namedtuple('Dimensions', ['dr', 'Nr', 'dtheta', 'Ntheta', 'Nframes'])
SysInfo = namedtuple('SysInfo', ['NL', 'NB', 'NBperTail', 'BoxArea', 'ExpBeadDensity', 'DrDtheta'])


class Symmetric_Site:
    """
    An aggregation of multiple binding sites on/in a protein/inclusion. User \
    defines the base_site Site object and then provides it to the \
    Symmetric_Site constructor.

    Attributes
    ----------
    temp : float
        The temperature of your system in K.

    Calculated Properties
    ---------------------
    symmetry : int
        The N-fold symmetry desired. I.E. 5 would yield 5 Sites.
    site_list : list
        The list of constituent Site objects that make up this Symmetric_Site.
    site_counts_histogram : numpy ndarray
        One-dimensional ndarray where the histogrammed ligand bead counts are \
        stored. e.g. [12, 5, 0, 0, 1, 0] would correspond to 12 frames having \
        zero beads in the Site, 5 frames having one bead in the Site, 0 frames \
        having 2, 3, or 5 beads in the site, and 1 frame having 4 beads in the \
        Site.
    bulk_counts_histogram : numpy ndarray
        One-dimensional ndarray where the histogrammed ligand bead counts are \
        stored. e.g. [12, 5, 0, 0, 1, 0] would correspond to 12 frames having \
        zero beads in the bulk patch, 5 frames having one bead in the patch, 0 \
        frames having 2, 3, or 5 beads in the patch, and 1 frame having 4 beads\
        in the patch.
    n_peak : int
        The mode of the bulk histogram. Indicates the cut-off for P_unocc.
    dG : float
        The binding affinity of the lipid for the Site, in kcal/mol.
    """

    def __init__(self, symmetry, base_site, counts_data):
        """
        Create a Symmetric_Site by creating clones of the base_site and \
        rotating them symmetrically around the origin.

        Parameters
        ----------
        symmetry : int
            The N-fold symmetry desired. I.E. 5 would yield 5 Sites.
        base_site : Site
            The original Site object that should be cloned and rotated.
        counts_data : ndarray
            The output from polarDensityBin, after having been read in and \
            processed by parse_tcl_data.

        Returns
        -------
        Symmetric_Site
            The Symmetric_Site you just created.
        list
            A list of the constituent Site objects that make up this \
            Symmetric_Site.

        """
        assert isinstance(symmetry, int), "symmetry must be an integer."
        assert isinstance(counts_data, np.ndarray), "counts_data must be a numpy ndarray."
        assert len(counts_data.shape) == 3, "counts_data must be a 3d array."
        Ntheta = counts_data.shape[2]
        assert Ntheta % symmetry == 0, f"The symmetry you specified does not \
            evenly distribute across the lattice you provided. Rerun with a \
            lattice that has theta bins divisible by {symmetry}. Current # of \
            theta bins is {Ntheta}."
        assert isinstance(base_site, Site), "base_site must be a Site."
        assert base_site.bin_coords is not None, "The base_site needs to be fully defined before creating a Symmetric_Site."
        self._symmetry = symmetry
        self._site_list = self._make_symmetric_sites(self, base_site, counts_data)
        assert len(self.site_list) == symmetry, "Number of Sites does not match symmetry."
        self.temp = base_site.temp

    @property
    def symmetry(self):
        """
        Tell me the symmetry, but don't let me change the symmetry.

        Returns
        -------
        int
            The N-fold symmetry of the Symmetric_Site.

        """
        return self._symmetry

    @property
    def site_list(self):
        """
        Tell me the site_list, but don't let me change the site_list.

        Returns
        -------
        list
            List of constituent Site objects that comprise this Symmetric_Site.

        """
        return self._site_list

    @property
    def site_counts_histogram(self):
        """
        Tell me the current counts, in histogram form, for the Symmetric_Site.

        Returns
        -------
        site_counts_histogram : numpy ndarray
            One-dimensional ndarray where the histogrammed ligand bead counts \
            are stored. e.g. [12, 5, 0, 0, 1, 0] would correspond to 12 frames \
            having zero beads in the Site, 5 frames having one bead in the \
            Site, 0 frames having 2, 3, or 5 beads in the site, and 1 frame \
            having 4 beads in the Site.

        """
        return _aggregate_site_counts_histograms(self.site_list)

    @property
    def bulk_counts_histogram(self):
        """
        Tell me the current counts, in histogram form, for the Symmetric_Site. \
        In practice, this is just the bulk_counts_histogram for the base_site.

        Returns
        -------
        bulk_counts_histogram : numpy ndarray
            One-dimensional ndarray where the histogrammed ligand bead counts \
            are stored. e.g. [12, 5, 0, 0, 1, 0] would correspond to 12 frames \
            having zero beads in the bulk patch, 5 frames having one bead in \
            the patch, 0 frames having 2, 3, or 5 beads in the patch, and 1 \
            frame having 4 beads in the patch.

        """
        return _check_bulk_counts_histogram(self.site_list)

    @property
    def n_peak(self):
        """
        Tell me what the n_peak is.

        Returns
        -------
        int
            The mode of the bulk distribution in a patch of membrane that has \
            equal accessible area to the site.

        """
        return calculate_hist_mode(self.bulk_counts_histogram)

    @property
    def dG(self):
        """
        Calculate the binding affinity of the lipid for this Symmetric_Site, \
        including the bulk correction factor dG_ref.

        Returns
        -------
        float
            The total binding affinity, in kcal/mol.

        """
        n_peak = self.n_peak
        dG_site = _calculate_dG(self.site_counts_histogram, n_peak, self.temp)
        dG_ref = _calculate_dG(self.bulk_counts_histogram, n_peak, self.temp)
        return dG_site - dG_ref

    def _make_symmetric_sites(self, base_site, counts_data):
        """
        Create identical sites to the base_site, rotated symmetrically around \
        the origin.

        Parameters
        ----------
        base_site : Site
            The Site object that you want to replicate symmetrically.
        counts_data : ndarray
            The output from polarDensityBin, after having been read in and \
            processed by parse_tcl_data.

        Returns
        -------
        site_list : list of Sites
            The list of all Sites that comprise this Symmetric_Site.

        """
        site_list = [base_site]
        Ntheta = counts_data.shape[2]
        base_site.name = base_site.name + '_1'
        bulk_counts_data = base_site.bulk_counts_histogram
        for site_number in range(1, self.symmetry):
            site_name = base_site.name + '_' + str(site_number + 1)
            new_site = Site(site_name, base_site.ligand, base_site.leaflet, base_site.temp)
            new_site.bin_coords = self._rotate_bin_coords(base_site.bin_coords, Ntheta, site_number)
            new_site.update_counts_histogram(False, counts_data)
            new_site.update_counts_histogram(True, bulk_counts_data)
            site_list.append(new_site)
        return site_list

    def _rotate_bin_coords(self, bin_coords, Ntheta, site_number):
        """
        Rotate the provided bin_coords around the circle.

        Parameters
        ----------
        bin_coords : list of tuples
            The bins that belong to this site in (r, theta) format. e.g. \
            [(2, 10), (2, 11), (2, 12)] would correspond to the 11th, 12th, and \
            13th theta bins in the 3rd radial bin from the origin. Bin coordinates \
            are zero-indexed by convention.
        Ntheta : int
            The number of theta bins in the circle.
        site_number : int
            Which constituent site this is.

        Returns
        -------
        rotated_bin_coords : list of tuples
            Should match input bin_coords in length and first tuple component; \
            second tuple components should all be shifted (rotated).

        """
        rotated_bin_coords = []
        for each_bin in bin_coords:
            r_bin, theta_bin = each_bin
            shift = Ntheta // self.symmetry
            rotated_theta_bin = theta_bin + shift * site_number
            if rotated_theta_bin >= Ntheta:
                rotated_theta_bin -= Ntheta
            rotated_bin_coords.append((r_bin, rotated_theta_bin))
        return rotated_bin_coords


class Site:
    """
    The basic class for a binding site on/in a protein/inclusion. User defines \
    Site with bin coordinates. Multiple symmetric Sites can be combined with \
    Symmetric_Site class. One Site (or Symmetric_Site) across multiple replicas \
    can be combined in Site_Across_Replicas class.

    Attributes
    ----------
    name : str
        The name of this Site. e.g. "Binding site 1," "Left anterior cleft," or \
        something else descriptive.
    ligand : str
        The name of the ligand that putatively binds this Site.
    leaflet : int
        If 1, outer leaflet. If 2, inner leaflet.
    temp : float
        The temperature of your system in K.

    Settable Properties
    -------------------
    bin_coords : list of tuples
        The bins that belong to this site in (r, theta) format. e.g. \
        [(2, 10), (2, 11), (2, 12)] would correspond to the 11th, 12th, and \
        13th theta bins in the 3rd radial bin from the origin. Bin coordinates \
        are zero-indexed by convention.

    Calculated Properties
    ---------------------
    site_counts_histogram : numpy ndarray
        One-dimensional ndarray where the histogrammed ligand bead counts are \
        stored. e.g. [12, 5, 0, 0, 1, 0] would correspond to 12 frames having \
        zero beads in the Site, 5 frames having one bead in the Site, 0 frames \
        having 2, 3, or 5 beads in the site, and 1 frame having 4 beads in the \
        Site.
    bulk_counts_histogram : numpy ndarray
        One-dimensional ndarray where the histogrammed ligand bead counts are \
        stored. e.g. [12, 5, 0, 0, 1, 0] would correspond to 12 frames having \
        zero beads in the bulk patch, 5 frames having one bead in the patch, 0 \
        frames having 2, 3, or 5 beads in the patch, and 1 frame having 4 beads\
        in the patch.
    n_peak : int
        The mode of the bulk histogram. Indicates the cut-off for P_unocc.
    dG : float
        The binding affinity of the lipid for the Site, in kcal/mol.
    """

    def __init__(self, name, ligand, leaflet, temp):
        """
        Create a Site object.

        Parameters
        ----------
        name : str
            The name of this Site. e.g. "Binding site 1," \
            "Left anterior cleft," or something else descriptive.
        ligand : str
            The name of the ligand that putatively binds this Site.
        leaflet : int
            If 1, outer leaflet. If 2, inner leaflet.
        temp : float
            The temperature of your system in K.
        """
        self.name = name
        self.ligand = ligand
        assert leaflet in [1, 2], "leaflet must be 1 or 2 (1 for outer leaflet or 2 for inner leaflet)"
        self.leaflet = leaflet
        self.temp = temp
        self._bin_coords = None
        self._site_counts_histogram = None
        self._bulk_counts_histogram = None

    @property
    def bin_coords(self):
        """
        Tell me what the bin_coords are. This is a getter function.

        Returns
        -------
        list of tuples
            The bins that belong to this site in (r, theta) format. e.g. \
            [(2, 10), (2, 11), (2, 12)] would correspond to the 11th, 12th, and \
            13th theta bins (starting at theta=0) in the 3rd radial bin from \
            the origin. Bin coordinates are zero-indexed by convention.

        """
        return self._bin_coords

    @bin_coords.setter
    def bin_coords(self, bin_coords):
        """
        Set bin_coords for this Site.

        Parameters
        ----------
        bin_coords : list of tuples
            The bins that belong to this site in (r, theta) format. e.g. \
            [(2, 10), (2, 11), (2, 12)] would correspond to the 11th, 12th, and \
            13th theta bins (starting at theta=0) in the 3rd radial bin from \
            the origin. Bin coordinates are zero-indexed by convention.

        Returns
        -------
        None.

        """
        assert isinstance(bin_coords, list), "bin_coords must be provided as a list"
        assert len(bin_coords) > 0, "bin_coords must have at least one bin"
        for bin_pair in bin_coords:
            assert isinstance(bin_pair, tuple), "bin_coords must be provided as a list of 2-tuples"
            assert len(bin_pair) == 2, f"bin_coords contains an invalid coordinate pair: {bin_pair}"
        self._bin_coords = bin_coords

    @property
    def site_counts_histogram(self):
        """
        Tell me the current counts, in histogram form, for the Site.

        Returns
        -------
        site_counts_histogram : numpy ndarray
            One-dimensional ndarray where the histogrammed ligand bead counts \
            are stored. e.g. [12, 5, 0, 0, 1, 0] would correspond to 12 frames \
            having zero beads in the Site, 5 frames having one bead in the \
            Site, 0 frames having 2, 3, or 5 beads in the site, and 1 frame \
            having 4 beads in the Site.

        """
        return self._site_counts_histogram

    @property
    def bulk_counts_histogram(self):
        """
        Tell me the current counts, in histogram form, for the Site.

        Returns
        -------
        bulk_counts_histogram : numpy ndarray
            One-dimensional ndarray where the histogrammed ligand bead counts \
            are stored. e.g. [12, 5, 0, 0, 1, 0] would correspond to 12 frames \
            having zero beads in the bulk patch, 5 frames having one bead in \
            the patch, 0 frames having 2, 3, or 5 beads in the patch, and 1 \
            frame having 4 beads in the patch.

        """
        return self._bulk_counts_histogram

    @property
    def n_peak(self):
        """
        Tell me what the n_peak is.

        Returns
        -------
        int
            The mode of the bulk distribution in a patch of membrane that has \
            equal accessible area to the site.

        """
        if self._bulk_counts_histogram is None:
            raise Exception("You need to update the bulk counts histogram first.")
        return calculate_hist_mode(self.bulk_counts_histogram)

    @property
    def dG(self):
        """
        Calculate the binding affinity of the lipid for this Site, including \
        the bulk correction factor dG_ref.

        Returns
        -------
        float
            The total binding affinity, in kcal/mol.

        """
        assert self.site_counts_histogram is not None, "You need to update the\
            site counts histogram first."
        assert self.bulk_counts_histogram is not None, "You need to add bulk \
            counts via update_counts_histogram(bulk=True, counts_data)."
        n_peak = self.n_peak
        assert n_peak is not None, "n_peak is missing."
        dG_site = _calculate_dG(self.site_counts_histogram, n_peak, self.temp)
        dG_ref = _calculate_dG(self.bulk_counts_histogram, n_peak, self.temp)
        return dG_site - dG_ref

    def _calculate_hist_mean(self, bulk=False):
        """
        Calculate the mean of the counts histogram.

        Parameters
        ----------
        bulk : boolean, optional
            If True, calculate the mean of the bulk distribution. Otherwise, \
            calculate the site mean. The default is False.

        Returns
        -------
        mean : float
            The mean of the counts histogram.

        """
        if bulk:
            hist = self.bulk_counts_histogram
        else:
            hist = self.site_counts_histogram
        total_N = np.sum(hist)
        sum_i = 0
        for i in range(len(hist)):
            sum_i += i * hist[i]
        mean = sum_i / total_N
        return mean

    def update_counts_histogram(self, bulk, counts_data):
        """
        Assign ligand bead counts to Site attribute "counts_histogram".

        Parameters
        ----------
        bulk : boolean
            If True, update the counts histogram for the bulk patch. If False,\
            update the counts histogram for the site.
        counts_data : ndarray
            If bulk=True, provide 1D nddarray containing bulk counts. \
            If bulk=False, provide the 3D ndarray containing binned counts.

        Returns
        -------
        None.

        """
        assert isinstance(counts_data, np.ndarray), "ndarray not supplied"
        if bulk:
            assert len(counts_data.shape) == 1, "Bulk counts data is not in the right format: {data}"
            bulk_hist = np.bincount(counts_data)
            self._bulk_counts_histogram = bulk_hist
            self._n_peak = calculate_hist_mode(self._bulk_counts_histogram)
        else:
            assert len(counts_data.shape) == 3, "Counts data is not in the right format: {data}"
            site_counts = self._fetch_site_counts(counts_data)
            site_hist = np.bincount(site_counts)
            self._site_counts_histogram = site_hist

    def _fetch_site_counts(self, binned_counts):
        """
        Create a 2D array where each row is a different bin within the site \
        and each column is a frame in the trajectory. Then sum over all the \
        rows to get site counts rather than bin counts.

        Parameters
        ----------
        binned_counts : ndarray
            3D ndarray containing the binned counts from polarDensityBin.tcl

        Returns
        -------
        site_counts : ndarray
            1D ndarray containing the site counts for each frame of the \
            trajectory.

        """
        stack = np.zeros(binned_counts.shape[0])
        for bin_tuple in self.bin_coords:
            r_bin, theta_bin = bin_tuple
            stack = np.vstack((stack, binned_counts[:, r_bin, theta_bin]))
        site_counts = np.sum(stack, axis=0)
        return site_counts

    def calculate_geometric_area(self, dr, dtheta):
        """
        Calculate the geometric area of the site.

        Parameters
        ----------
        dr : float
            The bin length in the radial dimension (Angstroms).
        dtheta : float
            The bin length in the azimuthal dimension (degrees).

        Returns
        -------
        float
            The geometric area of the site in square Angstroms.

        """
        area = 0
        for bin_tuple in self.bin_coords:
            area += _calculate_bin_area(bin_tuple[0], dr, dtheta)
        return area

    def predict_accessible_area(self, bulk_area, mode=False):
        """
        Predict the accessible area of the site. A reasonable method is to \
        multiply the area of the bulk patch you just analyzed by the ratio of\
        the means (or modes) for the site distribution and the bulk \
        distribution. This will put you in the ballpark of the bulk patch area.

        Parameters
        ----------
        bulk_area : float
            The area of the bulk patch previously analyzed in square Angstroms.
        mode : boolean
            If True, use the site and bulk modes rather than the means. Default\
            is False (use means).

        Returns
        -------
        predicted_accessible_area : float
            The area of the bulk patch you should analyze next to try and more\
            closely match the site distribution. Units are square Angstroms.

        """
        if self.ligand != "DPPC":
            warnings.append("Warning: This feature was originally built for use in a 100% DPPC test system. Proceed with caution.")
        if mode:
            site = calculate_hist_mode(self.site_counts_histogram)
            bulk = calculate_hist_mode(self.bulk_counts_histogram)
        else:
            site = self._calculate_hist_mean(bulk=False)
            bulk = self._calculate_hist_mean(bulk=True)
        predicted_accessible_area = bulk_area * (site / bulk)
        return predicted_accessible_area


def parse_tcl_dat_file(filepath, bulk):
    """
    Extract bin counts and bin dimensions from a TCL output .dat file.

    Parameters
    ----------
    filepath : Path
        The path to the TCL .dat output file.
    bulk : boolean
        If True, assume output is from do_get_counts. If False, assume \
        output is from polarDensityBin.

    Returns
    -------
    counts : ndarray
        If bulk is True, return a 1D ndarray. If bulk is False and avg is False,\
        return a 3D ndarray of integer counts. Dimension 0 is time, dimension 1\
        is r and dimension 2 is theta. IE the count for the 4th radial bin and \
        12th theta bin on the 33rd frame of the trajectory would be \
        counts[32, 3, 11].
    dimensions : namedtuple or None
        If bulk is True, return None. If bulk is False, return the bin \
        dimensions in r and theta as well as the number of frames inside a \
        namedtuple.
    system_info : namedtuple or None
        The information taken from the header row of the .dat file specified in\
        filepath. If bulk is True, returns None.

    """
    if isinstance(filepath, str):
        filepath = Path(filepath)
    elif not isinstance(filepath, Path):
        raise Exception("Must provide the path to the .dat file.")
    assert filepath.exists(), f"Could not find {filepath}"
    assert filepath.is_file(), f"This is not recognized as a file {filepath}"
    assert filepath.suffixes[-1] == '.dat', "You must provide the .dat file output from VMD."
    if bulk:
        return np.loadtxt(filepath).astype(int).flatten(), None, None
    else:
        unrolled_data = np.loadtxt(filepath, skiprows=1)
        system_info = _parse_system_info(np.loadtxt(filepath, comments=None, max_rows=1, delimiter=',', dtype=str))
        grid_dims = _calculate_grid_dimensions(unrolled_data)
        counts = _package_counts(unrolled_data, grid_dims).squeeze()
        return counts, grid_dims, system_info


def calculate_density(avg_counts, grid_dims):
    """
    Calculate the average bead density in each bin.

    Parameters
    ----------
    avg_counts : ndarray
        2D array of average bead counts per bin in simulation. [r,theta] format.
    grid_dims : namedtuple
        The Dimensions object corresponding to this system.

    Returns
    -------
    density : ndarray
        The average bead density in each bin.

    """
    assert isinstance(grid_dims, tuple), "grid_dims must be a Dimensions namedtuple."
    assert isinstance(avg_counts, np.ndarray), "avg_counts must be an ndarray"
    assert len(avg_counts.shape) == 2, "avg_counts must be a 2D array."
    area = _calculate_lattice_areas(grid_dims)
    density = avg_counts / area
    return density


def _calculate_grid_dimensions(unrolled_data):
    """
    Return the length and number of bins in each dimension, as well as the \
    number of frames in the .dat file.

    Parameters
    ----------
    unrolled_data : ndarray
        The data, as it is read in directly from the .dat file output by \
        polarDensityBin.

    Returns
    -------
    grid_dims : namedtuple
        Contains dr, number of r bins, dtheta, number of theta bins, and number\
        of frames contained in file.

    """
    dr = unrolled_data[0, 1] - unrolled_data[0, 0]
    dthetadeg = unrolled_data[0, 2]
    dtheta = dthetadeg * (np.pi / 180.0)
    nframes = _calculate_nframes(unrolled_data[:, 0])
    Ntheta = int(round(360 / dthetadeg))
    assert Ntheta == unrolled_data.shape[1] - 4, f"Something went wrong with the theta dimensions parser. dtheta={dtheta}, Ntheta={Ntheta}"
    Nr = len(unrolled_data[:, 0]) / nframes
    assert Nr - int(Nr) == 0, f"Something went wrong with the r dimensions parser. dr={dr}, Nr={Nr}"
    Nr = int(Nr)
    grid_dims = Dimensions(dr, Nr, dtheta, Ntheta, nframes)
    return grid_dims


def _package_counts(unrolled_data, grid_dims):
    """
    Package the unrolled data into a 3d array in [time, r, theta] format.

    Parameters
    ----------
    unrolled_data : ndarray
        2d array output from polarDensityBin.
    grid_dims : namedtuple
        Contains bin dimensions and number of frames.

    Returns
    -------
    counts : ndarray
        3d array of counts in [time, r, theta] format.

    """
    nframes = grid_dims.Nframes

    # chop off the first few columns
    unrolled_counts = unrolled_data[:, 3:-1]

    # 'sideways' because it is in [r, time, theta] format at first
    sideways_counts = np.zeros((grid_dims.Nr, nframes, grid_dims.Ntheta))
    for i in range(grid_dims.Nr):
        sideways_counts[i, :, :] = unrolled_counts[(nframes * i):(nframes * (i + 1)), :]

    # swap axes to put it in [time, r, theta] format
    counts = np.swapaxes(sideways_counts, 0, 1)
    return counts


def _calculate_lattice_areas(grid_dims):
    """
    Calculate the area of each bin in a polar lattice.

    Parameters
    ----------
    grid_dims : namedtuple
        Dimensions namedtuple that corresponds to your system.

    Returns
    -------
    areas : ndarray
        2D array of bin areas in [r, theta] format.

    """
    areas = np.zeros((grid_dims.Nr, grid_dims.Ntheta))
    for radial_ring in range(grid_dims.Nr):
        areas[radial_ring, :] = _calculate_bin_area(radial_ring, grid_dims.dr, grid_dims.dtheta)
    return areas


def _calculate_bin_area(r_bin, dr, dtheta):
    """
    Calculate the area of the polar bin.

    Parameters
    ----------
    r_bin : int
        Which radial bin is this? Zero-indexed.
    dr : float
        The radial bin length in Angstroms.
    dtheta : float
        The azimuthal bin length in degrees.

    Returns
    -------
    area : float
        The bin area in square Angstroms.

    """
    bin_radial_midpoint = (r_bin * dr) + (0.5 * dr)
    area = dr * dtheta * bin_radial_midpoint
    return area


def _calculate_nframes(r_values):
    """
    Calculate how many frames are in a trajectory from the first column of the \
    .dat file produced by the TCL script.

    Parameters
    ----------
    r_values : ndarray
        The beginning r values for each bin in your system.

    Returns
    -------
    nframes : int
        The number of frames in your trajectory.

    """
    match_value = r_values[0]
    for i in range(len(r_values)):
        if r_values[i] != match_value:
            return i
    return len(r_values)


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
    NL, NB, BoxArea, ExpBeadDensity, NBperTail, DrDtheta = dat_file_header
    NL = _isolate_number_from_header_string(NL)
    NB = _isolate_number_from_header_string(NB)
    BoxArea = _isolate_number_from_header_string(BoxArea)
    ExpBeadDensity = _isolate_number_from_header_string(ExpBeadDensity)
    NBperTail = _isolate_number_from_header_string(NBperTail)
    DrDtheta = _isolate_number_from_header_string(DrDtheta)
    sysInfo = SysInfo(NL, NB, NBperTail, BoxArea, ExpBeadDensity, DrDtheta)
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
    else:
        return float(right_side.split()[0])


def _calculate_dG(counts_histogram, n_peak, temp):
    """
    Calculate the delta G for bulk or site.

    Parameters
    ----------
    counts_histogram : np.ndarray
        The site or bulk counts histogram attribute of a Site.
    n_peak : int
        The mode of the bulk counts histogram.
    temp : float
        The temperature of your system in K.

    Returns
    -------
    delta_G : float
        The binding affinity for the site or bulk, in kcal/mol.

    """
    minus_RT = -1.0 * temp * constants.R / 4184.  # 4184 converts J to kcal
    P_unnoc = _calculate_P_unnoc(counts_histogram, n_peak)
    delta_G = minus_RT * np.log((1 - P_unnoc) / P_unnoc)
    return delta_G


def _calculate_P_unnoc(counts_histogram, n_peak):
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
    P_unnoc = counts_histogram[:n_peak + 1] / total_N
    P_occ = counts_histogram[n_peak + 1:] / total_N
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
        warnings.append(f"Warning: More than one peak identified ({mode}), using first peak ({mode[0]})")
        mode = mode[0]
    if mode == 0:
        if nonzero:
            warnings.append("Warning: found an experimental mode of 0 for site, using second highest peak")
            mode = np.argmax(histogram[1:]) + 1
    return mode


def _aggregate_site_counts_histograms(site_list):
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
    first_site = site_list.pop(0)
    counts = first_site.site_counts_histogram
    for site in site_list:
        counts_to_add = site.site_counts_histogram
        # make sure the two ndarrays have same shape; resize the smaller one if not
        if counts.shape[0] > counts_to_add.shape[0]:
            counts_to_add = counts_to_add.resize(counts.shape)
        elif counts.shape[0] < counts_to_add.shape[0]:
            counts = counts.resize(counts_to_add.shape)
        counts += counts_to_add
    return counts


def _check_bulk_counts_histogram(site_list):
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
    first_site = site_list.pop(0)
    bulk = first_site.bulk_counts_histogram
    for site in site_list:
        assert bulk.all() == site.bulk_counts_histogram.all(), "One or more sites have different bulk histograms. This shouldn't be possible."
    return bulk


def outline_site(ax, site, grid_dims):
    """
    Draw an outline around each bin in this Site.

    Parameters
    ----------
    ax : matplotlib.pyplot Axes object
        The Axes object you want to draw this on.
    site : Site or Symmetric_Site
        The site you want to outline.
    grid_dims : namedtuple
        Contains dr, number of r bins, dtheta, number of theta bins, and number\
        of frames contained in file.

    Returns
    -------
    ax : matplotlib.pyplot Axes object
        The Axes object you want to draw this on.

    """
    if isinstance(site, Site):
        for each_bin in site.bin_coords:
            ax = outline_bin(ax, each_bin, grid_dims)
    elif isinstance(site, Symmetric_Site):
        for each_site in site.site_list:
            for each_bin in each_site.bin_coords:
                ax = outline_bin(ax, each_bin, grid_dims)
    else:
        Exception("site must be a Site or Symmetric_Site.")
    return ax


def outline_bin(ax, bin_coords, grid_dims):
    """
    Draw an outline around this bin.

    Parameters
    ----------
    ax : matplotlib.pyplot Axes object
        The Axes object you want to draw this on.
    bin_coords : tuple
        The tuple of bin coordinates stored in (r_bin, theta_bin) format.
    grid_dims : namedtuple
        Contains dr, number of r bins, dtheta, number of theta bins, and number\
        of frames contained in file.

    Returns
    -------
    ax : matplotlib.pyplot Axes object
        The Axes object you want to draw this on.

    """
    assert isinstance(bin_coords, tuple), "bin_coords must be a tuple of bin coordinates."
    dr, _, dtheta, _, _ = grid_dims
    start_theta = bin_coords[1] * dtheta
    end_theta = start_theta + dtheta
    inner_r = bin_coords[0] * dr
    outer_r = inner_r + dr
    theta_range = np.linspace(start_theta, end_theta, 100)
    ax.fill_between(theta_range, inner_r, outer_r, facecolor=(0, 0, 0, 0), edgecolor='k')
    return ax


def make_custom_colormap():
    """
    Make a custom colormap for plotting.

    Returns
    -------
    my_cmap : matplotlib colormap
        The colormap to use when plotting.

    """
    depleted = plt.colormaps['RdGy_r']
    enriched = plt.colormaps['bwr']
    newcolors = np.concatenate([depleted(np.linspace(0.35, 0.5, 128)), enriched(np.linspace(0.5, 1, 128))])
    my_cmap = ListedColormap(newcolors)
    return my_cmap


def create_heatmap_figure_and_axes(numlipids, cmap, v_vals, figwidth, figheight):
    """
    Create the heatmap figure and enough axes to accommodate all the lipids and\
    leaflets.

    Parameters
    ----------
    numlipids : int
        The number of lipids you intend to plot.
    cmap : matplotlib ListedColormap
        A custom colormap.
    v_vals : tuple
        The (vmin, vmid, and vmax)
    figwidth : int, optional
        Figure width.
    figheight : int, optional
        Figure height.

    Returns
    -------
    fig : matplotlib Fig object
        The figure you just created.
    matplotlib Axes objects
        The polar projection axes that were created.

    """
    vmin, vmid, vmax = v_vals
    fig = plt.figure(figsize=(figwidth, figheight))
    gs = gridspec.GridSpec(numlipids, 2, figure=fig, wspace=0.15, hspace=0.15)
    for gridbox in range(numlipids * 2):
        ax = plt.subplot(gs[gridbox], projection='polar')
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.21, 1, 0.5, 0.02])
    sm = mpl.cm.ScalarMappable(cmap=cmap)
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
    cbar.set_ticks(np.linspace(0, 1, 5))
    cbar.ax.set_xticklabels([vmin, (vmin + vmid) / 2, vmid, (vmid + vmax) / 2, vmax])
    return fig, fig.axes


def bin_prep(bin_info):
    """
    Configure the arrays needed for plotting polar heatmaps.

    Parameters
    ----------
    bin_info : namedtuple
        Contains Nr, Ntheta, dr, and dtheta information.

    Returns
    -------
    list
        The two numpy ndarrays needed for plotting a heatmap.

    """
    r_vals = np.linspace(0, bin_info.Nr * bin_info.dr, bin_info.Nr + 1)
    theta_vals = np.linspace(0, 2 * np.pi, bin_info.Ntheta + 1)
    r_vals, theta_vals = np.meshgrid(r_vals, theta_vals, indexing='ij')
    return [r_vals, theta_vals]


def plot_heatmap(ax, data, grid, cmap, v_vals):
    vmin, vmid, vmax = v_vals
    norm = MidpointNormalize(midpoint=vmid, vmin=vmin, vmax=vmax)
    ax.grid(False)
    plt.axis('off')
    radius, theta = grid
    ax.pcolormesh(theta, radius, data, cmap=cmap, norm=norm, zorder=0, edgecolors='face', linewidth=0)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    return ax


# This class comes from Liam Sharp and could potentially be rewritten to be
# more clear.
class MidpointNormalize(Normalize):
    """
    Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value).

    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """

    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
