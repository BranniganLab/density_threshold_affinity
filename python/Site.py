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
        return self.calculate_hist_mode(bulk=True, nonzero=False)

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
        return self.calculate_dG(bulk=False) - self.calculate_dG(bulk=True)

    def calculate_dG(self, bulk=False):
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
        P_unnoc = self.calculate_P_unnoc(bulk=bulk)
        delta_G = minus_RT * np.log((1 - P_unnoc) / P_unnoc)
        return delta_G

    def calculate_hist_mode(self, bulk=False, nonzero=False):
        """
        Calculate the mode of the counts histogram.

        Parameters
        ----------
        nonzero : boolean
            If True, in the case of mode=0, use the second highest peak instead.\
            This is necessary when estimating the predicted accessible area.

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

    def calculate_hist_mean(self, bulk=False):
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

    def update_counts_histogram(self, filepath, bulk=False):
        """
        Assign ligand bead counts to Site attribute "counts_histogram".

        Parameters
        ----------
        filepath : Path
            The path to the .dat file to be parsed.
        bulk : boolean
            If True, update the counts histogram for the bulk patch. If False,\
            update the counts histogram for the site. Default is False.

        Returns
        -------
        None.

        """
        if bulk:
            counts = np.loadtxt(filepath).astype(int).flatten()
            hist = np.bincount(counts)
            self._bulk_counts_histogram = hist
        else:
            counts = self.fetch_site_counts(filepath)
            hist = np.bincount(counts)
            self._site_counts_histogram = hist

    def calculate_P_unnoc(self, bulk=False):
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


    def fetch_site_counts(self, filepath):
        """
        

        Parameters
        ----------
        filepath : TYPE
            DESCRIPTION.

        Returns
        -------
        counts : TYPE
            DESCRIPTION.

        """
        return counts

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

    def predict_accessible_area(self, bulk_area):
        """
        

        Parameters
        ----------
        bulk_area : TYPE
            DESCRIPTION.

        Returns
        -------
        predicted_accessible_area : TYPE
            DESCRIPTION.

        """
        site_mean = self.calculate_hist_mean(bulk=False)
        bulk_mean = self.calculate_hist_mean(bulk=True)
        predicted_accessible_area = bulk_area * (site_mean / bulk_mean)
        return predicted_accessible_area
