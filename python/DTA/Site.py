#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 13:51:00 2024.

@author: js2746
"""
import numpy as np
from DTA.utils import calculate_hist_mode, calculate_hist_mean, calculate_dG
from DTA.density import calculate_bin_area


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
    leaflet_id : int
        If 1, outer leaflet. If 2, inner leaflet.
    temperature : float
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
    site_counts_over_time : numpy ndarray
        One-dimensional ndarray containing the total count for each frame.
    n_peak : int
        The mode of the bulk histogram. Indicates the cut-off for P_unocc.
    dG : float
        The binding affinity of the lipid for the Site, in kcal/mol.
    """

    def __init__(self, name, leaflet_id, temperature):
        """
        Create a Site object.

        Parameters
        ----------
        name : str
            The name of this Site. e.g. "Binding site 1," \
            "Left anterior cleft," or something else descriptive.
        leaflet_id : int
            If 1, outer leaflet. If 2, inner leaflet.
        temperature : float
            The temperature of your system in K.
        """
        self.name = name
        assert leaflet_id in [1, 2], "leaflet_id must be 1 or 2 (1 for outer leaflet or 2 for inner leaflet)"
        self.leaflet_id = leaflet_id
        self.temperature = temperature
        self._bin_coords = None
        self._site_counts_histogram = None
        self._bulk_counts_histogram = None
        self._site_counts_over_time = None

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
    def site_counts_over_time(self):
        """
        Tell me the counts on each frame of the trajectory for this Site.

        Returns
        -------
        counts_over_time : numpy ndarray
            One-dimensional ndarray containing the total count for each frame.

        """
        return self._site_counts_over_time

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
        dG_site = calculate_dG(self.site_counts_histogram, n_peak, self.temperature)
        dG_ref = calculate_dG(self.bulk_counts_histogram, n_peak, self.temperature)
        return dG_site - dG_ref

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
            assert len(counts_data.shape) == 1, f"Bulk counts data is not in the right format: {counts_data}"
            bulk_hist = np.bincount(counts_data)
            self._bulk_counts_histogram = bulk_hist
            self._n_peak = calculate_hist_mode(self._bulk_counts_histogram)
        else:
            assert len(counts_data.shape) == 3, f"Counts data is not in the right format: {counts_data}"
            counts_data = counts_data.astype(int)
            self._site_counts_over_time = self._fetch_site_counts(counts_data)
            self._site_counts_histogram = np.bincount(self._site_counts_over_time)

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
            area += calculate_bin_area(bin_tuple[0], dr, dtheta)
        return area

    def predict_accessible_area(self, bulk_area, mode=True):
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
            is True. If False, use means (untested feature).

        Returns
        -------
        predicted_accessible_area : float
            The area of the bulk patch you should analyze next to try and more\
            closely match the site distribution. Units are square Angstroms.

        """
        if mode:
            site = calculate_hist_mode(self.site_counts_histogram)
            bulk = calculate_hist_mode(self.bulk_counts_histogram)
        else:
            site = calculate_hist_mean(self.site_counts_histogram)
            bulk = calculate_hist_mean(self.bulk_counts_histogram)
        predicted_accessible_area = bulk_area * (site / bulk)
        return predicted_accessible_area

    def _fetch_site_counts(self, binned_counts):
        """
        Create a 2D array where each row is a different bin within the site \
        and each column is a frame in the trajectory. Then sum over all the \
        rows to get site counts over time.

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
        return site_counts.astype(int)
