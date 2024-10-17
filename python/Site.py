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

Dimensions = namedtuple('Dimensions', ['dr', 'Nr', 'dtheta', 'Ntheta', 'Nframes'])


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
        13th theta bins (starting at theta=0) in the 3rd radial bin from the \
        origin. Bin coordinates are zero-indexed by convention.

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
        self._n_peak = None
        self._dG = None
        return self

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
    def bin_coords(self, bin_coords, filepath=None):
        """
        Set bin_coords and in the process recalculate the counts histogram.

        Parameters
        ----------
        bin_coords : list of tuples
            The bins that belong to this site in (r, theta) format. e.g. \
            [(2, 10), (2, 11), (2, 12)] would correspond to the 11th, 12th, and \
            13th theta bins (starting at theta=0) in the 3rd radial bin from \
            the origin. Bin coordinates are zero-indexed by convention.
        filepath : str or Path or None
            The path to the .dat file output by PolarDensity_for_DTA.tcl. If \
            provided, the site counts histogram will automatically update. \
            Default value is None.

        Returns
        -------
        None.

        """
        assert isinstance(bin_coords, list), "bin_coords must be provided as a list"
        assert len(self.bin_coords) > 0, "bin_coords must have at least one bin"
        for bin_pair in bin_coords:
            assert isinstance(bin_pair, tuple), "bin_coords must be provided as a list of 2-tuples"
            assert len(bin_pair) == 2, f"bin_coords contains an invalid coordinate pair: {bin_pair}"
        self._bin_coords = bin_coords
        if filepath is not None:
            if isinstance(filepath, str):
                filepath = Path(filepath)
            self._site_counts_histogram = self.update_counts_histogram(filepath, bulk=False)

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
        return self._calculate_hist_mode(bulk=True, nonzero=False)

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
        return self._calculate_dG(bulk=False) - self._calculate_dG(bulk=True)

    def _calculate_dG(self, bulk=False):
        """
        Calculate the delta G. If bulk is True, this value is dG_ref. If bulk \
        is False, this value is dG_site.

        Parameters
        ----------
        bulk : boolean
            Is this the bulk patch? The default is False.

        Returns
        -------
        delta_G : float
            The binding affinity for the site or bulk, in kcal/mol.

        """
        minus_RT = -1.0 * self.temp * constants.R / 4184.  # kcal/mol
        P_unnoc = self._calculate_P_unnoc(bulk=bulk)
        delta_G = minus_RT * np.log((1 - P_unnoc) / P_unnoc)
        return delta_G

    def _calculate_hist_mode(self, bulk=False, nonzero=False):
        """
        Calculate the mode of the counts histogram.

        Parameters
        ----------
        nonzero : boolean
            If True, in the case of mode=0, use the second highest peak instead.\
            This may be necessary when estimating the predicted accessible area.

        Returns
        -------
        mode : int
            The mode of the distribution.

        """
        if bulk:
            hist = self.bulk_counts_histogram
        else:
            hist = self.site_counts_histogram
        mode = np.argmax(hist)
        if len(np.shape(mode)) != 0:
            warnings.append(f"Warning: More than one peak identified ({mode}), using first peak ({mode[0]})")
            mode = mode[0]
        if mode == 0:
            if nonzero:
                warnings.append(f"Warning: found an experimental mode of 0 for site '{self.name}', using second highest peak")
                mode = np.argmax(hist[1:]) + 1
        return mode

    def _calculate_hist_mean(self, bulk=False):
        """
        xtz.

        Parameters
        ----------
        bulk : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        mean : TYPE
            DESCRIPTION.

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

    def update_counts_histogram(self, bulk, data):
        """
        Assign ligand bead counts to Site attribute "counts_histogram".

        Parameters
        ----------
        bulk : boolean
            If True, update the counts histogram for the bulk patch. If False,\
            update the counts histogram for the site.
        data : ndarray or str or Path
            If bulk=True, provide the path (str or Path) to the .dat file from \
            do_get_counts.tcl. If bulk=False, provide the 3D ndarray containing\
            binned counts information.

        Returns
        -------
        None.

        """
        if bulk:
            if isinstance(data, str):
                filepath = Path(data)
            elif isinstance(data, Path):
                filepath = data
            else:
                raise Exception("Must provide the path to the bulk counts.")
            assert filepath.exists(), f"Could not find {filepath}"
            assert filepath.is_file(), f"This is not recognized as a file {filepath}"
            assert filepath.suffixes[1] == '.dat', "You must provide the .dat file output from do_get_counts.tcl in VMD."
            bulk_counts = np.loadtxt(filepath).astype(int).flatten()
            bulk_hist = np.bincount(bulk_counts)
            self._bulk_counts_histogram = bulk_hist
        else:
            assert isinstance(data, np.ndarray), "ndarray not supplied"
            assert len(data.shape) == 3, "This ndarray is not 3 dimensional"
            site_counts = self._fetch_site_counts(data)
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

    def _calculate_P_unnoc(self, bulk=False):
        """
        Calculate the probability that the site is unoccupied by ligand. This \
        is defined by the portion of the histogram <= n_peak divided by the \
        total histogram.

        Parameters
        ----------
        bulk : boolean
            If True, calculate P_unnoc,bulk. Else, calculate P_unnoc,site.

        Returns
        -------
        P_unnoc : float
            The probability that the site is unoccupied by ligand.

        """
        if bulk:
            counts = self.bulk_counts_histogram
        else:
            counts = self.site_counts_histogram
        total_N = np.sum(counts)
        P_unnoc = counts[:self.n_peak + 1] / total_N
        P_occ = counts[self.n_peak + 1:] / total_N
        assert math.isclose(P_occ + P_unnoc, 1, abs_tol=0.01), f"Probabilities do not sum to one. Current sum: {P_unnoc + P_occ}"
        return P_unnoc

    def calculate_geometric_area(self, dr, dtheta):
        """
        Calculate the geometric area of the site.

        Parameters
        ----------
        dr : float
            The bin length in the radial dimension.
        dtheta : float
            The bin length in the azimuthal dimension.

        Returns
        -------
        float
            The geometric area of the site.

        """
        area = 0
        for bin_tuple in self.bin_coords:
            bin_radial_midpoint = (bin_tuple[0] * dr) + (0.5 * dr)
            area += dr * dtheta * bin_radial_midpoint
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
            The area of the bulk patch previously analyzed.
        mode : boolean
            If True, use the site and bulk modes rather than the means. Default\
            is False (use means).

        Returns
        -------
        predicted_accessible_area : float
            The area of the bulk patch you should analyze next to try and more\
            closely match the site distribution.

        """
        if self.ligand != "DPPC":
            warnings.append("Warning: This feature was originally built for use in a 100% DPPC test system. Proceed with caution.")
        if mode:
            site = self._calculate_hist_mode(bulk=False)
            bulk = self._calculate_hist_mode(bulk=True)
        else:
            site = self._calculate_hist_mean(bulk=False)
            bulk = self._calculate_hist_mean(bulk=True)
        predicted_accessible_area = bulk_area * (site / bulk)
        return predicted_accessible_area


def calculate_nframes(r_values):
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
    return None


def parse_tcl_dat_trajectory(filepath, bulk):
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
        If bulk is True, return a 1D ndarray. If bulk is False, return a 3D \
        ndarray of integer counts. Dimension 0 is time, dimension 1 is r and \
        dimension 2 is theta. IE the count for the 4th radial bin and 12th \
        theta bin on the 33rd frame of the trajectory would be counts[32, 3, 11].
    dimensions : namedtuple or None
        If bulk is True, return None. If bulk is False, return the bin \
        dimensions in r and theta as well as the number of frames inside a \
        namedtuple.

    """
    if bulk:
        return np.loadtxt(filepath).astype(int).flatten(), None
    else:
        unrolled_data = np.loadtxt(filepath, skiprows=1, delimiter=' ')

        # calculate polar lattice dimensions, nframes
        dr = unrolled_data[0, 1] - unrolled_data[0, 0]
        dtheta = unrolled_data[0, 2]
        nframes = calculate_nframes(unrolled_data[:, 0])
        Ntheta = int(round(360 / dtheta))
        assert Ntheta == unrolled_data.shape[1] - 3, f"Something went wrong with the theta dimensions parser. dtheta={dtheta}, Ntheta={Ntheta}"
        Nr = len(unrolled_data[:, 0]) / nframes
        assert Nr - int(Nr) == 0, f"Something went wrong with the r dimensions parser. dr={dr}, Nr={Nr}"
        Nr = int(Nr)
        grid_dims = Dimensions(dr, Nr, dtheta, Ntheta, nframes)

        # package bin counts into 3d ndarray in [time, r, theta] format
        unrolled_counts = unrolled_data[:, 3:].astype(int)
        sideways_counts = np.zeros((Nr, nframes, Ntheta))
        for i in range(Nr):
            sideways_counts[i, :, :] = unrolled_counts[(nframes * i):(nframes * (i + 1)), :]
        counts = np.swapaxes(sideways_counts, 0, 1)
        return counts, grid_dims
