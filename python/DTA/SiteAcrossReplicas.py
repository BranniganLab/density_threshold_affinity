#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 17:20:36 2025.

@author: js2746
"""
import numpy as np
from DTA.Site import Site
from DTA.SymmetricSite import SymmetricSite
from DTA.utils import calculate_hist_mode, calculate_hist_mean, calculate_dG, aggregate_site_counts_histograms, check_bulk_counts_histogram


class SiteAcrossReplicas:
    """
    An aggregation of a single binding site across multiple replicas. User \
    defines the base_site SymmetricSite and/or Site object first (including \
    setting the bin_coords!) and then provides it to the Site_Across_Replicas \
    constructor.

    Attributes
    ----------
    temperature : float
        The temperature of your system in K.

    Calculated Properties
    ---------------------
    name : str
        The name of the Site. Will be inherited from base_site.
    get_site_list : list
        The list of constituent Site or SymmetricSite objects that make up \
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
        SymmetricSites and/or Sites that comprise this SiteAcrossReplicas.
    """

    def __init__(self, replica_list, base_site):
        """
        Create a Site_Across_Replicas.

        Parameters
        ----------
        replica_list : list
            List of all the different replicas you want to include.
        base_site : Site or SymmetricSite
            The original Site object that should be used across replicas. Can \
            be Site or SymmetricSite.

        """
        if isinstance(base_site, Site):
            assert base_site.bin_coords is not None, "The base_site needs to be fully defined before creating a Site_Across_Replicas."
        elif not isinstance(base_site, SymmetricSite):
            raise ValueError("base_site must be a Site or SymmetricSite")
        self.name = base_site.name
        self._site_list = self._make_sites_across_replicas(base_site, replica_list)
        assert len(self.get_site_list) == len(replica_list), "Number of Sites does not match number of replicas."
        self.temperature = base_site.temperature

    def __iter__(self):
        """Iterate through the site_list."""
        yield from self.get_site_list

    @property
    def get_site_list(self):
        """
        Tell me the site_list, but don't let me change the site_list.

        Returns
        -------
        list
            List of constituent Site objects that comprise this SymmetricSite.

        """
        return self._site_list

    @property
    def site_counts_histogram(self):
        """
        Tell me the current counts, in histogram form, for the SymmetricSite.

        Returns
        -------
        site_counts_histogram : numpy ndarray
            One-dimensional ndarray where the histogrammed ligand bead counts \
            are stored. e.g. [12, 5, 0, 0, 1, 0] would correspond to 12 frames \
            having zero beads in the Site, 5 frames having one bead in the \
            Site, 0 frames having 2, 3, or 5 beads in the site, and 1 frame \
            having 4 beads in the Site.

        """
        return aggregate_site_counts_histograms(self.get_site_list)

    @property
    def bulk_counts_histogram(self):
        """
        Tell me the current counts, in histogram form, for the SymmetricSite. \
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
        return check_bulk_counts_histogram(self.get_site_list)

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
        Calculate the binding affinity of the lipid for this SymmetricSite, \
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
        constituent Sites that comprise this SymmetricSite.

        Returns
        -------
        float
            The standard deviation

        """
        dGs = []
        for site in self.get_site_list:
            dGs.append(site.dG)
        return np.std(np.array(dGs))

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
        for site in self.get_site_list:
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

    def _make_sites_across_replicas(self, base_site, replica_list):
        """
        Create identical sites to the base_site, across multiple replicas.

        Parameters
        ----------
        base_site : Site
            The Site object that you want to replicate symmetrically.
        replica_list : list
            The list of replicas.

        Returns
        -------
        site_list : list of Sites
            The list of all Sites that comprise this SymmetricSite.

        """
        assert isinstance(replica_list, list), "replica_list must be a list"
        name = base_site.name
        site_list = []
        for site_number, _ in enumerate(replica_list):
            site_name = name + '_rep' + str(site_number + 1)
            if isinstance(base_site, Site):
                new_site = Site(site_name, base_site.leaflet_id, base_site.temperature)
                new_site.bin_coords = base_site.bin_coords
            elif isinstance(base_site, SymmetricSite):
                new_single_site = Site(site_name, base_site.get_site_list[0].leaflet_id, base_site.get_site_list[0].temperature)
                new_single_site.bin_coords = base_site.get_site_list[0].bin_coords
                new_site = SymmetricSite(base_site.symmetry, new_single_site, base_site._Ntheta)
            new_site.update_counts_histogram(bulk=False, counts_data=replica_list[site_number])
            site_list.append(new_site)
        return site_list
