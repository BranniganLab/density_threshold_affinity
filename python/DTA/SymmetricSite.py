#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 16:53:18 2024.

@author: js2746
"""
import numpy as np
from DTA.Site import Site
from DTA.utils import calculate_hist_mode, calculate_hist_mean, calculate_dG, aggregate_site_counts_histograms, check_bulk_counts_histogram


class SymmetricSite:
    """
    An aggregation of multiple binding sites on/in a protein/inclusion. User \
    defines the base_site Site object first (including setting the bin_coords!)\
    and then provides it to the SymmetricSite constructor.

    Attributes
    ----------
    temperature : float
        The temperature of your system in K.

    Calculated Properties
    ---------------------
    name : str
        The name of the Site. Will be inherited from base_site.
    symmetry : int
        The N-fold symmetry desired. I.E. 5 would yield 5 Sites.
    ntheta_in_lattice : int
        The number of azimuthal bins in the lattice (not just how many are in \
        the site!)
    bin_coords : list of tuples
        The bins that belong to this site in (r, theta) format. e.g. \
        [(2, 10), (2, 11), (2, 12)] would correspond to the 11th, 12th, and \
        13th theta bins in the 3rd radial bin from the origin. Bin coordinates \
        are zero-indexed by convention.
    get_site_list : list
        The list of constituent Site objects that make up this SymmetricSite.
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
    site_counts_over_time : list
        First item is a one-dimensional ndarray containing the total count for \
        each frame. Second item is list of one-dimensional ndarrays containing \
        total counts for each Site that comprises this SymmetricSite.
    n_peak : int
        The mode of the bulk histogram. Indicates the cut-off for P_unocc.
    dG : float
        The binding affinity of the lipid for the Site, in kcal/mol.
    dG_std : float
        The standard deviation of the mean binding affinity for the Sites that\
        comprise this SymmetricSite.
    """

    def __init__(self, symmetry, base_site, Ntheta):
        """
        Create a SymmetricSite by creating clones of the base_site and \
        rotating them symmetrically around the origin.

        Parameters
        ----------
        symmetry : int
            The N-fold symmetry desired. I.E. 5 would yield 5 Sites.
        base_site : Site
            The original Site object that should be cloned and rotated.

        """
        assert isinstance(symmetry, int), "symmetry must be an integer."
        assert Ntheta % symmetry == 0, "This symmetry does not evenly divide \
            across the number of theta bins."
        assert isinstance(base_site, Site), "base_site must be a Site."
        assert base_site.bin_coords is not None, "The base_site needs to be fully defined before creating a SymmetricSite."
        self.name = base_site.name
        self._symmetry = symmetry
        self._ntheta_in_lattice = Ntheta
        self._site_list = self._make_symmetric_sites(base_site, Ntheta)
        assert len(self.get_site_list) == symmetry, "Number of Sites does not match symmetry."
        self.temperature = base_site.temperature

    def __iter__(self):
        """Iterate through the site_list."""
        yield from self.get_site_list

    @property
    def ntheta_in_lattice(self):
        """Tell me how many azimuthal bins there are in the lattice (not this \
        site)."""
        return self._ntheta_in_lattice

    @property
    def symmetry(self):
        """
        Tell me the symmetry, but don't let me change the symmetry.

        Returns
        -------
        int
            The N-fold symmetry of the SymmetricSite.

        """
        return self._symmetry

    @property
    def bin_coords(self):
        """
        Generate one list of bin coordinate tuples corresponding to all the \
        bins inside this SymmetricSite. Necessary for outline_site.

        Returns
        -------
        bin_coords_list : list of tuples
            The bins that belong to this SymmetricSite in (r, theta) format. \
            e.g. [(2, 10), (2, 11), (2, 12)] would correspond to the 11th, \
            12th, and 13th theta bins (starting at theta=0) in the 3rd radial \
            bin from the origin. Bin coordinates are zero-indexed by convention.

        """
        bin_coords_list = []
        for site in self.get_site_list:
            site_coords = site.bin_coords
            for each_bin in site_coords:
                bin_coords_list.append(each_bin)
        return bin_coords_list

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
    def site_counts_over_time(self):
        """
        Tell me the counts on each frame of the trajectory for this SymmetricSite\
        and its constituent Sites.

        Returns
        -------
        agg_counts_over_time : numpy ndarray
            One-dimensional ndarray containing the aggregated total \
            SymmetricSite counts for each frame.
        ind_counts_over_time : list of numpy ndarrays
            List containing the site_counts_over_time for each individual Site \
            within this SymmetricSite.

        """
        ind_counts_over_time = []
        for site in self.get_site_list:
            ind_counts_over_time.append(site.site_counts_over_time)
        agg_counts_over_time = np.sum(np.array(ind_counts_over_time), axis=0)
        return agg_counts_over_time, ind_counts_over_time

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
            The list of all Sites that comprise this SymmetricSite.

        """
        name = base_site.name
        base_site.name = name + '_1'
        site_list = [base_site]
        for site_number in range(1, self.symmetry):
            site_name = name + '_' + str(site_number + 1)
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
