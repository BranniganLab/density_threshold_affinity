#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 13:51:00 2024.

@author: js2746
"""
from __future__ import annotations
from typing import Literal
import numpy as np
from dta.utils import calculate_hist_mode, calculate_hist_mean, calculate_dG
from dta.bin_logic import PolarBinGrid, BinAddress


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
    grid : PolarBinGrid
        Contains lattice information.
    leaflet_id : int
        If 1, outer leaflet. If 2, inner leaflet.
    temperature : float
        The temperature of your system in K.

    Settable Properties
    -------------------
    bin_coords : frozenset of BinAddress
        The bins that belong to this site in (r, theta) format. e.g. \
        {(2, 10), (2, 11), (2, 12)} would correspond to the 11th, 12th, and \
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

    def __init__(self, name: str, grid: PolarBinGrid, leaflet_id: Literal[1, 2], temperature: float):
        """
        Create a Site object.

        Parameters
        ----------
        name : str
            The name of this Site. e.g. "Binding site 1," \
            "Left anterior cleft," or something else descriptive.
        grid : PolarBinGrid
            Contains lattice information.
        leaflet_id : int
            If 1, outer leaflet. If 2, inner leaflet.
        temperature : float
            The temperature of your system in K.
        """
        self.name = name
        if not isinstance(grid, PolarBinGrid):
            raise TypeError("grid must be a PolarBinGrid object.")
        self.grid = grid
        if leaflet_id not in [1, 2]:
            raise ValueError("leaflet_id must be 1 or 2 (1 for outer leaflet or 2 for inner leaflet)")
        if not np.isfinite(temperature):
            raise ValueError("temperature must be finite.")
        if temperature <= 0:
            raise ValueError("temperature must be positive.")
        self.leaflet_id = leaflet_id
        self.temperature = temperature
        self._bin_coords = None
        self._site_counts_histogram = None
        self._bulk_counts_histogram = None

    @property
    def bin_coords(self) -> frozenset[BinAddress] | None:
        """
        Return the bins defining this Site.

        The returned collection is immutable. Assign a new collection to
        ``bin_coords`` to change the Site definition.

        Returns
        -------
        frozenset of BinAddress
            The bins that belong to this site in (r, theta) format. e.g. \
            [(2, 10), (2, 11), (2, 12)] would correspond to the 11th, 12th, and \
            13th theta bins (starting at theta=0) in the 3rd radial bin from \
            the origin. Bin coordinates are zero-indexed by convention.

        """
        return self._bin_coords

    @bin_coords.setter
    def bin_coords(
        self,
        bin_addresses: (
            list[BinAddress]
            | tuple[BinAddress, ...]
            | set[BinAddress]
            | frozenset[BinAddress]
        )
    ) -> None:
        """
        Validate and replace the bins defining this Site.

        Parameters
        ----------
        bin_addresses : list, tuple, set, or frozenset of BinAddress
            The bins that belong to this site in (r, theta) format. e.g.
            [(2, 10), (2, 11), (2, 12)] would correspond to the 11th, 12th, and
            13th theta bins (starting at theta=0) in the 3rd radial bin from
            the origin. Bin coordinates are zero-indexed by convention.

        Returns
        -------
        None.

        """
        if not isinstance(bin_addresses, (list, tuple, set, frozenset)):
            raise TypeError("bin_addresses must be provided as a list, tuple, set, or frozenset")
        if not bin_addresses:
            raise ValueError("bin_addresses cannot be empty.")
        validated = set()
        for item in bin_addresses:
            address = item if isinstance(item, BinAddress) else BinAddress(*item)
            if not 0 <= address.r_index < self.grid.r.n_bins:
                raise IndexError(
                    f"Radial bin {address.r_index} out of range."
                )
            if not 0 <= address.theta_index < self.grid.theta.n_bins:
                raise IndexError(
                    f"Angular bin {address.theta_index} out of range."
                )
            validated.add(address)

        new_bin_coords = frozenset(validated)

        if new_bin_coords != self._bin_coords:
            self._bin_coords = new_bin_coords
            self._site_counts_histogram = None
            self._bulk_counts_histogram = None

    @property
    def site_counts_histogram(self) -> np.ndarray:
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
    def bulk_counts_histogram(self) -> np.ndarray:
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
    def n_peak(self) -> float:
        """
        Tell me what the n_peak is.

        Returns
        -------
        int
            The mode of the bulk distribution in a patch of membrane that has \
            equal accessible area to the site.

        """
        if self._bulk_counts_histogram is None:
            raise AttributeError("You need to update the bulk counts histogram first.")
        return calculate_hist_mode(self.bulk_counts_histogram)

    @property
    def dG(self) -> float:
        """
        Calculate the binding affinity of the lipid for this Site, including \
        the bulk correction factor dG_ref.

        Returns
        -------
        float
            The total binding affinity, in kcal/mol.

        """
        if self.site_counts_histogram is None:
            raise RuntimeError("You need to update the site counts histogram before calculating dG.")
        if self.bulk_counts_histogram is None:
            raise RuntimeError("You need to update the bulk counts histogram before calculating dG.")
        n_peak = self.n_peak
        if n_peak is None:
            raise RuntimeError("n_peak is missing.")
        dG_site = calculate_dG(self.site_counts_histogram, n_peak, self.temperature)
        dG_ref = calculate_dG(self.bulk_counts_histogram, n_peak, self.temperature)
        return dG_site - dG_ref

    def copy(self, new_name: str) -> Site:
        """Return a copy of this site with empty histograms."""
        copy = Site(new_name, self.grid, self.leaflet_id, self.temperature)
        if self.bin_coords is not None:
            copy.bin_coords = self.bin_coords
        return copy

    def update_site_counts_histogram(self, counts_data: np.ndarray) -> None:
        """
        Assign ligand bead counts to Site attribute "site_counts_histogram".

        Parameters
        ----------
        counts_data : np.ndarray
            The 3D ndarray containing binned site counts.

        Returns
        -------
        None.

        """
        if self.bin_coords is None:
            raise RuntimeError(
                "Cannot update the site counts histogram before bin_coords "
                "have been defined."
            )
        counts_data = self._validate_counts_data(counts_data=counts_data, expected_ndim=3)
        if counts_data.shape[-2:] != (self.grid.r.n_bins, self.grid.theta.n_bins):
            raise ValueError(f"""
            counts_data is the wrong shape for this lattice.
            {counts_data.shape} != {(self.grid.r.n_bins, self.grid.theta.n_bins)}
            """)
        site_counts = self._fetch_site_counts(counts_data)
        site_hist = np.bincount(site_counts)
        self._site_counts_histogram = site_hist

    def update_bulk_counts_histogram(self, counts_data: np.ndarray) -> None:
        """
        Assign ligand bead counts to Site attribute "bulk_counts_histogram".

        Parameters
        ----------
        counts_data : ndarray
            The 1D nddarray containing bulk counts.

        Returns
        -------
        None.

        """
        if self.bin_coords is None:
            raise RuntimeError(
                "Cannot update the bulk counts histogram before bin_coords "
                "have been defined."
            )
        counts_data = self._validate_counts_data(counts_data=counts_data, expected_ndim=1)
        bulk_hist = np.bincount(counts_data)
        self._bulk_counts_histogram = bulk_hist

    @staticmethod
    def _validate_counts_data(
        counts_data: np.ndarray,
        expected_ndim: int,
    ) -> np.ndarray:
        """
        Validate and normalize frame-resolved molecular count data.

        Confirm that ``counts_data`` is a NumPy array with the dimensionality
        required by the calling histogram-update method. The array must have a
        real numeric dtype and contain only finite, integer-valued,
        non-negative values. Integer-valued floating-point arrays, such as
        ``[0.0, 1.0, 2.0]``, are accepted because some input parsers represent
        discrete counts as floats. Fractional values are rejected rather than
        silently truncated.

        After validation, return a new array containing the same values cast
        to NumPy's default integer dtype. The returned array is suitable for
        indexing operations and for use with ``numpy.bincount``; the input
        array is not modified.

        Parameters
        ----------
        counts_data : numpy.ndarray
            Frame-resolved count data to validate. Bulk counts are expected to
            be one-dimensional, while spatially binned site counts are
            expected to be three-dimensional.
        expected_ndim : int
            Number of dimensions required by the calling method.

        Returns
        -------
        numpy.ndarray
            A newly allocated integer array with the same shape and count
            values as ``counts_data``.

        Raises
        ------
        TypeError
            If ``counts_data`` is not a NumPy array or does not have a real
            numeric dtype.
        ValueError
            If the array has the wrong number of dimensions or contains NaN,
            infinity, fractional values, or negative counts.
        """
        if not isinstance(counts_data, np.ndarray):
            raise TypeError(
                "counts_data must be provided as a NumPy ndarray."
            )

        if counts_data.ndim != expected_ndim:
            raise ValueError(
                f"counts_data must be {expected_ndim}D; "
                f"received shape {counts_data.shape}."
            )

        if counts_data.shape[0] == 0:
            raise ValueError("counts_data must include at least one frame.")

        if (
            not np.issubdtype(counts_data.dtype, np.number)
            or np.issubdtype(counts_data.dtype, np.complexfloating)
        ):
            raise TypeError(
                "counts_data must contain real numeric values; "
                f"received dtype {counts_data.dtype}."
            )

        if not np.all(np.isfinite(counts_data)):
            raise ValueError(
                "counts_data must contain only finite values."
            )

        if not np.all(counts_data == np.floor(counts_data)):
            raise ValueError(
                "counts_data must contain integer-valued counts."
            )

        if np.any(counts_data < 0):
            raise ValueError(
                "counts_data cannot contain negative counts."
            )

        max_supported_count = np.iinfo(np.intp).max

        if np.any(counts_data > max_supported_count):
            raise ValueError(
                "counts_data contains values too large for this platform; "
                f"the maximum supported count is {max_supported_count}."
            )

        return counts_data.astype(np.intp)

    def calculate_geometric_area(self) -> float:
        """
        Calculate the geometric area of the site.

        Returns
        -------
        float
            The geometric area of the site in square Angstroms.

        """
        area = 0.0
        for bin_address in self.bin_coords:
            area += self.grid.calc_bin_area(bin_address)
        return area

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

    def _fetch_site_counts(self, binned_counts: np.ndarray) -> np.ndarray[int]:
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
