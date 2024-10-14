#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 13:44:55 2024.

@author: js2746
"""

import numpy as np
import warnings
from utils import fetch_bin_count, fetch_bin_area
from nougat import Membrane


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
    bin_coords : list of tuples
        The bins that belong to this site in (r, theta) format. e.g. \
        [(2, 10), (2, 11), (2, 12)] would correspond to the 11th, 12th, and \
        13th theta bins (starting at theta=0) in the 3rd radial bin from the \
        origin. Bin coordinates are zero-indexed by convention.
    counts_hist : numpy ndarray
        One-dimensional ndarray where the histogrammed ligand bead counts are \
        stored. e.g. [12, 5, 0, 0, 1, 0] would correspond to 12 frames having \
        zero beads in the Site, 5 frames having one bead in the Site, 0 frames \
        having 2, 3, or 5 beads in the site, and 1 frame having 4 beads in the \
        Site.
    n_peak : int
        The mode of the bulk histogram. Indicates the cut-off for P_unocc.
    P_unocc : float
        The probability that the Site is unoccupied by ligand.
    dG : float
        The binding affinity of ligand for this Site, in kcal/mol.
    """

    def __init__(self, name, ligand, leaflet):
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
        """
        self.name = name
        self.ligand = ligand
        assert leaflet in [1, 2], "leaflet must be 1 or 2 (1 for outer leaflet or 2 for inner leaflet)"
        self.leaflet = leaflet
        self._bin_coords = None
        self._counts_hist = None
        self._n_peak = None
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
    def bin_coords(self, bin_coords, membrane_obj):
        """
        Set bin_coords and in the process recalculate the counts histogram.

        Parameters
        ----------
        bin_coords : list of tuples
            The bins that belong to this site in (r, theta) format. e.g. \
            [(2, 10), (2, 11), (2, 12)] would correspond to the 11th, 12th, and \
            13th theta bins (starting at theta=0) in the 3rd radial bin from \
            the origin. Bin coordinates are zero-indexed by convention.
        membrane_obj : Membrane
            A nougat Membrane class object that contains the counts for the \
            given ligand and leaflet.

        Returns
        -------
        None.

        """
        assert isinstance(bin_coords, list), "bin_coords must be provided as a list"
        assert len(self.bin_coords) > 0, "bin_coords must have at least one bin"
        assert isinstance(bin_coords[0], tuple), "bin_coords must be provided as a list of 2-tuples"
        for bin_pair in bin_coords:
            assert len(bin_pair) == 2, f"bin_coords contains an invalid coordinate pair: {bin_pair}"
        assert isinstance(membrane_obj, Membrane), "membrane_obj must be a valid nougat Membrane object"
        self._bin_coords = bin_coords
        self._counts_hist = self.update_counts_hist(membrane_obj)

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
        return self._n_peak

    @n_peak.setter
    def n_peak(self, n_peak):
        """
        Set n_peak and in the process recalculate the P_unnoc and dG.

        Parameters
        ----------
        n_peak : int
            The mode of the bulk distribution in a patch of membrane that has \
            equal accessible area to the site.

        Returns
        -------
        None.

        """
        if n_peak is not None:
            assert isinstance(n_peak, int), "Non-integer n_peak is not supported."
            # assert self._bin_coords is not None, "Please define the site (bin) coordinates before adding n_peak."
        self._n_peak = n_peak

    @property
    def site_mode(self):
        """
        Calculate the mode of the distribution. If the mode is zero, use the\
        next-closest value.

        Returns
        -------
        site_mode : int
            The non-zero mode of the site distribution.

        """
        site_mode = np.argmax(self._counts_hist)
        if site_mode == 0:
            warnings.append(f"Warning: found an experimental mode of 0 for site '{self.name}', using second highest peak")
            site_mode = np.argmax(self._counts_hist[1:]) + 1
        return site_mode

    def update_counts_hist(self, membrane_obj):
        """
        Assign ligand bead counts to Site attribute "counts_hist".

        Parameters
        ----------
        None.

        Returns
        -------
        None.

        """
        count = np.zeros(membrane_obj.grid_dims['Nframes'])
        for bin_tuple in self.bin_coords:
            count += fetch_bin_count(bin_tuple)
        self._counts_hist = count

    def calculate_geometric_area(self, membrane_obj):
        """
        Calculate the geometric area of the site.

        Parameters
        ----------
        membrane_obj : nougat Membrane
            A nougat Membrane class object that contains the counts for the \
            given ligand and leaflet.

        Returns
        -------
        float
            The geometric area of the site.

        """
        grid_dims_dict = membrane_obj.grid_dims
        dr = grid_dims_dict['d1']
        dtheta = grid_dims_dict['d2']
        area = 0
        for bin_tuple in self.bin_coords:
            bin_radial_midpoint = (bin_tuple[0] * dr) + (0.5 * dr)
            area += dr * dtheta * bin_radial_midpoint
        return area
