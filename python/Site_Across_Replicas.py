#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 16:53:18 2024.

@author: js2746
"""
from Site import Site
from Symmetric_Site import Symmetric_Site
import numpy as np
from utils import calculate_hist_mode, calculate_hist_mean, calculate_dG


class Site_Across_Replicas:
    """
    An aggregation of a single binding site across multiple replicas. User \
    defines the base_site Symmetric_Site and/or Site object first (including \
    setting the bin_coords!) and then provides it to the Site_Across_Replicas \
    constructor.

    Attributes
    ----------
    temperature : float
        The temperature of your system in K.

    Calculated Properties
    ---------------------
    site_list : list
        The list of constituent Site or Symmetric_Site objects that make up \
        this Site_Across_Replicas.
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
    dG_std : float
        The standard deviation of the mean binding affinity for the \
        Symmetric_Sites and/or Sites that comprise this Site_Across_Replicas.
    """

    def __init__(self, replica_list, base_site):
        """
        Create a Site_Across_Replicas.

        Parameters
        ----------
        replica_list : list
            List of all the different replicas you want to include.
        base_site : Site or Symmetric_Site
            The original Site object that should be used across replicas. Can \
            be Site or Symmetric_Site.

        """
        if isinstance(base_site, Site):
            assert base_site.bin_coords is not None, "The base_site needs to be fully defined before creating a Site_Across_Replicas."
        elif not isinstance(base_site, Symmetric_Site):
            raise Exception("base_site must be a Site or Symmetric_Site")
        self._site_list = self._make_sites_across_replicas(base_site, replica_list)
        assert len(self.site_list) == len(replica_list), "Number of Sites does not match number of replicas."
        self.temperature = base_site.temperature

    def __iter__(self):
        """Iterate through the site_list."""
        for site in self.site_list:
            yield site

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
        dG_site = calculate_dG(self.site_counts_histogram, n_peak, self.temperature)
        dG_ref = calculate_dG(self.bulk_counts_histogram, n_peak, self.temperature)
        return dG_site - dG_ref

    @property
    def dG_std(self):
        """
        Calculate the standard deviation of the delta G values across the \
        constituent Sites that comprise this Symmetric_Site.

        Returns
        -------
        float
            The standard deviation

        """
        dGs = []
        for site in self.site_list:
            dGs.append(site.dG)
        return np.std(dGs)

    def update_counts_histogram(self, bulk, counts_data):
        """
        Update the counts histograms for all constituent Sites.

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
        for site in self.site_list:
            site.update_counts_histogram(bulk, counts_data)

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
            is True. If False use means instead of modes (untested feature).

        Returns
        -------
        predicted_accessible_area : float
            The area of the bulk patch you should analyze next to try and more\
            closely match the site distribution. Units are in square Angstroms.

        """
        if mode:
            site = calculate_hist_mode(self.site_counts_histogram)
            bulk = calculate_hist_mode(self.bulk_counts_histogram)
        else:
            site = calculate_hist_mean(self.site_counts_histogram)
            bulk = calculate_hist_mean(self.bulk_counts_histogram)
        predicted_accessible_area = bulk_area * (site / bulk)
        return predicted_accessible_area

    def _make_symmetric_sites(self, base_site, Ntheta):
        """
        Create identical sites to the base_site, rotated symmetrically around \
        the origin.

        Parameters
        ----------
        base_site : Site
            The Site object that you want to replicate symmetrically.
        Ntheta : int
            The total number of theta bins in the circle.

        Returns
        -------
        site_list : list of Sites
            The list of all Sites that comprise this Symmetric_Site.

        """
        base_site.name = base_site.name + '_1'
        site_list = [base_site]
        for site_number in range(1, self.symmetry):
            site_name = base_site.name + '_' + str(site_number + 1)
            new_site = Site(site_name, base_site.leaflet_id, base_site.temperature)
            new_site.bin_coords = self._rotate_bin_coords(base_site.bin_coords, Ntheta, site_number)
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
    hist_lengths = []
    for site in site_list:
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
    first_site = site_list[0]
    bulk = first_site.bulk_counts_histogram.copy()
    for site in site_list[1:]:
        assert bulk.all() == site.bulk_counts_histogram.all(), "One or more sites have different bulk histograms. This shouldn't be possible."
    return bulk
