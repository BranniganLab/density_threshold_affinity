#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 16:53:18 2024.

@author: js2746
"""
from __future__ import annotations
import numpy as np
from dta.Site import Site
from dta.bin_logic import BinAddress
from dta.utils import calculate_hist_mode, calculate_hist_mean, calculate_dG, aggregate_site_counts_histograms, check_bulk_counts_histogram


class SymmetricSite:
    """
    An aggregation of multiple binding sites on/in a protein/inclusion. User \
    defines the base_site Site object first (including setting the bin_coords!)\
    and then provides it to the SymmetricSite constructor.

    Attributes
    ----------
    name: str
        A user-generated name for this site. Inherits from base_site.
    temperature : float
        The temperature of your system in K. Inherits from base_site.
    grid : PolarBinGrid
        Contains lattice information. Inherits from base_site.

    Calculated Properties
    ---------------------
    symmetry : int
        The N-fold symmetry desired. I.E. 5 would yield 5 Sites.
    bin_coords : set of BinAddress
        The bins that belong to this site as BinAddress objects. Bin indices
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
    n_peak : int
        The mode of the bulk histogram. Indicates the cut-off for P_unocc.
    dG : float
        The binding affinity of the lipid for the Site, in kcal/mol.
    dG_std : float
        The standard deviation of the mean binding affinity for the Sites that\
        comprise this SymmetricSite.
    """

    def __init__(self, symmetry: int, base_site: Site):
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
        if isinstance(symmetry, bool) or not isinstance(symmetry, int):
            raise TypeError("symmetry must be an integer.")
        if symmetry < 2:
            raise ValueError("symmetry must be greater than one.")
        if not isinstance(base_site, Site):
            raise TypeError("base_site must be a Site object.")
        if base_site.grid.theta.n_bins % symmetry != 0:
            raise ValueError("This polar lattice is not evenly divisible by the symmetry number provided.")
        if base_site.bin_coords is None:
            raise ValueError("The base_site needs to be fully defined before creating a SymmetricSite.")
        self.name = base_site.name
        self._symmetry = symmetry
        self._site_list = self._make_symmetric_sites(base_site)
        assert len(self.get_site_list) == symmetry, "Number of Sites does not match symmetry."
        self.temperature = base_site.temperature
        self.grid = base_site.grid

    def __iter__(self):
        """Iterate through the site_list."""
        yield from self.get_site_list

    @property
    def symmetry(self) -> int:
        """
        Tell me the symmetry, but don't let me change the symmetry.

        Returns
        -------
        int
            The N-fold symmetry of the SymmetricSite.

        """
        return self._symmetry

    @property
    def bin_coords(self) -> set[BinAddress]:
        """
        Generate one list of BinAddress[es] corresponding to all the \
        bins inside this SymmetricSite. Necessary for outline_site.

        Returns
        -------
        set of BinAddress
            The bins that belong to this SymmetricSite as BinAddress objects.
            Bin indices are zero-indexed by convention.

        """
        bin_coords_list = []
        for site in self.get_site_list:
            site_coords = site.bin_coords
            for each_bin in site_coords:
                bin_coords_list.append(each_bin)
        return set(bin_coords_list)

    @property
    def get_site_list(self) -> list[Site]:
        """
        Tell me the site_list, but don't let me change the site_list.

        Returns
        -------
        list
            List of constituent Site objects that comprise this SymmetricSite.

        """
        return self._site_list

    @property
    def site_counts_histogram(self) -> np.ndarray:
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
    def bulk_counts_histogram(self) -> np.ndarray:
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
    def n_peak(self) -> int:
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
    def dG(self) -> float:
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
    def dG_std(self) -> float:
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

    def copy(self, new_name: str) -> SymmetricSite:
        """Return a copy of this SymmetricSite with empty histograms."""
        base_site_copy = self.get_site_list[0].copy(new_name)
        symm_site_copy = SymmetricSite(self.symmetry, base_site_copy)
        return symm_site_copy

    def update_site_counts_histogram(self, counts_data: np.ndarray) -> None:
        """
        Update the site counts histograms for all constituent Sites.

        Parameters
        ----------
        counts_data : ndarray
            The 3D ndarray containing binned site counts.

        Returns
        -------
        None.

        """
        for site in self.get_site_list:
            site.update_site_counts_histogram(counts_data)

    def update_bulk_counts_histogram(self, counts_data: np.ndarray) -> None:
        """
        Update the bulk counts histograms for all constituent Sites.

        Parameters
        ----------
        counts_data : ndarray
            The 1D ndarray containing bulk counts.

        Returns
        -------
        None.

        """
        for site in self.get_site_list:
            site.update_bulk_counts_histogram(counts_data)

    def predict_accessible_area(self, bulk_area: float, mode: bool = True) -> float:
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

    def _make_symmetric_sites(self, base_site: Site) -> list[Site]:
        """
        Create identical sites to the base_site, rotated symmetrically around \
        the origin.

        Parameters
        ----------
        base_site : Site
            The Site object that you want to replicate symmetrically.

        Returns
        -------
        site_list : list of Sites
            The list of all Sites that comprise this SymmetricSite.

        """
        site_list = []
        for site_number in range(self.symmetry):
            site_name = f"{base_site.name}_{site_number + 1}"
            new_site = base_site.copy(site_name)
            new_site.bin_coords = self._rotate_bin_coords(
                base_site.bin_coords,
                base_site.grid.theta.n_bins,
                site_number,
            )
            site_list.append(new_site)
        self._check_for_overlapping_sites(site_list)
        return site_list

    def _rotate_bin_coords(self, bin_coords: list[BinAddress], n_theta: int, site_number: int) -> list[BinAddress]:
        """
        Rotate the provided bin_coords around the circle.

        Parameters
        ----------
        bin_coords : list of BinAddress
            The bins that belong to this site in (r, theta) format. e.g. \
            [(2, 10), (2, 11), (2, 12)] would correspond to the 11th, 12th, and \
            13th theta bins in the 3rd radial bin from the origin. Bin coordinates \
            are zero-indexed by convention.
        n_theta : int
            The number of angular bins in the lattice.
        site_number : int
            Which constituent site this is.

        Returns
        -------
        rotated_bin_coords : list of BinAddress
            Should match input bin_coords in length and radial bin; theta bins
            should all be shifted (rotated).

        """
        if site_number >= self.symmetry:
            raise ValueError(f"Site number {site_number} is greater than total symmetry ({self.symmetry})")
        rotated_bin_coords = []
        for each_bin in bin_coords:
            r_bin, theta_bin = each_bin
            shift = n_theta // self.symmetry
            rotated_theta_bin = theta_bin + shift * site_number
            if rotated_theta_bin >= n_theta:
                rotated_theta_bin -= n_theta
            rotated_bin_coords.append(BinAddress(r_bin, rotated_theta_bin))
        return rotated_bin_coords

    def _check_for_overlapping_sites(self, site_list: list[Site]) -> None:
        """Raise an error when two constituent Sites share one or more bins."""
        for first_index, first_site in enumerate(site_list):
            for second_site in site_list[first_index + 1:]:
                overlap = first_site.bin_coords & second_site.bin_coords

                if overlap:
                    raise ValueError(
                        f"{first_site.name} and {second_site.name} overlap in "
                        f"{len(overlap)} bin(s): {sorted(overlap)}"
                    )
