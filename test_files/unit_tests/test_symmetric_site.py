#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for the SymmetricSite class."""

import numpy as np
import pytest

from dta.Site import Site
from dta.SymmetricSite import SymmetricSite
from dta.bin_logic import BinAddress, PolarBinGrid
from dta.utils import calculate_dG


@pytest.fixture
def grid():
    """Return a grid divisible into four symmetric sectors."""
    return PolarBinGrid(r_min=0.0, r_max=2.0, n_r=2, n_theta=8)


@pytest.fixture
def base_site(grid):
    """Return a fully defined base Site."""
    site = Site("binding site", grid, leaflet_id=2, temperature=320)
    site.bin_coords = {BinAddress(0, 7), BinAddress(1, 0)}
    return site


@pytest.fixture
def symmetric_site(base_site):
    """Return a fourfold SymmetricSite."""
    return SymmetricSite(symmetry=4, base_site=base_site)


def _counts_for_symmetric_site():
    """Return counts that differ among the four constituent sites."""
    counts = np.zeros((5, 2, 8), dtype=int)
    counts[:, 0, 7] = [0, 1, 1, 2, 2]
    counts[:, 0, 1] = [0, 1, 1, 2, 2]
    counts[:, 0, 3] = [0, 1, 2, 2, 2]
    counts[:, 0, 5] = [1, 1, 1, 2, 3]
    return counts


def _populate_histograms(symmetric_site):
    """Populate all constituent Site histograms."""
    symmetric_site.update_counts_histogram(False, _counts_for_symmetric_site())
    symmetric_site.update_counts_histogram(True, np.array([0, 1, 1, 2, 2]))


def test_init_creates_named_rotated_sites_and_inherits_attributes(base_site):
    symmetric = SymmetricSite(symmetry=4, base_site=base_site)

    assert symmetric.name == "binding site"
    assert symmetric.symmetry == 4
    assert symmetric.grid is base_site.grid
    assert symmetric.temperature == base_site.temperature
    assert [site.name for site in symmetric] == [
        "binding site_1",
        "binding site_2",
        "binding site_3",
        "binding site_4",
    ]
    assert [site.bin_coords for site in symmetric] == [
        {BinAddress(0, 7), BinAddress(1, 0)},
        {BinAddress(0, 1), BinAddress(1, 2)},
        {BinAddress(0, 3), BinAddress(1, 4)},
        {BinAddress(0, 5), BinAddress(1, 6)},
    ]


def test_iter_yields_constituent_sites(symmetric_site):
    assert list(symmetric_site) == symmetric_site.get_site_list
    assert all(isinstance(site, Site) for site in symmetric_site)


def test_bin_coords_returns_union_of_constituent_bins(symmetric_site):
    assert symmetric_site.bin_coords == {
        BinAddress(0, 1),
        BinAddress(0, 3),
        BinAddress(0, 5),
        BinAddress(0, 7),
        BinAddress(1, 0),
        BinAddress(1, 2),
        BinAddress(1, 4),
        BinAddress(1, 6),
    }


def test_init_rejects_non_integer_symmetry(base_site):
    with pytest.raises(TypeError, match="symmetry must be an integer"):
        SymmetricSite(2.0, base_site)


@pytest.mark.parametrize("symmetry", [0, -1])
def test_init_rejects_nonpositive_symmetry(base_site, symmetry):
    with pytest.raises(ValueError, match="symmetry must be positive"):
        SymmetricSite(symmetry, base_site)


def test_init_rejects_non_site_base():
    with pytest.raises(TypeError, match="base_site must be a Site"):
        SymmetricSite(2, object())


def test_init_rejects_symmetry_that_does_not_divide_grid(grid):
    site = Site("site", grid, 1, 300)
    site.bin_coords = [(0, 0)]

    with pytest.raises(ValueError, match="not evenly divisible"):
        SymmetricSite(3, site)


def test_init_requires_defined_base_site(grid):
    site = Site("site", grid, 1, 300)

    with pytest.raises(ValueError, match="base_site needs to be fully defined"):
        SymmetricSite(2, site)


def test_update_counts_histogram_updates_every_constituent_site(symmetric_site):
    counts = _counts_for_symmetric_site()
    bulk_counts = np.array([0, 1, 1, 2, 2])

    symmetric_site.update_counts_histogram(False, counts)
    symmetric_site.update_counts_histogram(True, bulk_counts)

    for site in symmetric_site:
        assert site.site_counts_histogram is not None
        np.testing.assert_array_equal(site.bulk_counts_histogram, np.array([1, 2, 2]))


def test_site_counts_histogram_aggregates_constituent_histograms(symmetric_site):
    symmetric_site.update_counts_histogram(False, _counts_for_symmetric_site())

    expected = np.zeros(4)
    for site in symmetric_site:
        expected[:len(site.site_counts_histogram)] += site.site_counts_histogram
    np.testing.assert_array_equal(symmetric_site.site_counts_histogram, expected)


def test_site_counts_histogram_requires_populated_sites(symmetric_site):
    with pytest.raises(AssertionError, match="do not have counts associated"):
        _ = symmetric_site.site_counts_histogram


def test_bulk_counts_histogram_and_n_peak(symmetric_site):
    symmetric_site.update_counts_histogram(True, np.array([0, 1, 1, 2, 2]))

    np.testing.assert_array_equal(symmetric_site.bulk_counts_histogram, np.array([1, 2, 2]))
    assert symmetric_site.n_peak == 1


def test_dg_uses_aggregated_site_and_bulk_histograms(symmetric_site):
    _populate_histograms(symmetric_site)

    expected = calculate_dG(symmetric_site.site_counts_histogram, 1, 320)
    expected -= calculate_dG(symmetric_site.bulk_counts_histogram, 1, 320)
    assert symmetric_site.dG == pytest.approx(expected)


def test_dg_std_is_standard_deviation_across_constituent_sites(symmetric_site):
    _populate_histograms(symmetric_site)

    expected = np.std(np.array([site.dG for site in symmetric_site]))
    assert symmetric_site.dG_std == pytest.approx(expected)


def test_predict_accessible_area_uses_modes(symmetric_site):
    _populate_histograms(symmetric_site)

    # Both aggregated-site and bulk histograms have mode 1.
    assert symmetric_site.predict_accessible_area(100.0) == pytest.approx(100.0)


def test_predict_accessible_area_can_use_means(symmetric_site):
    _populate_histograms(symmetric_site)

    site_histogram = symmetric_site.site_counts_histogram
    site_mean = sum(i * count for i, count in enumerate(site_histogram)) / sum(site_histogram)
    bulk_histogram = symmetric_site.bulk_counts_histogram
    bulk_mean = sum(i * count for i, count in enumerate(bulk_histogram)) / sum(bulk_histogram)
    expected = 100.0 * site_mean / bulk_mean
    assert symmetric_site.predict_accessible_area(100.0, mode=False) == pytest.approx(expected)


def test_copy_preserves_definition_but_rebuilds_empty_constituent_sites(symmetric_site):
    _populate_histograms(symmetric_site)

    copied = symmetric_site.copy("copy")

    assert copied is not symmetric_site
    assert copied.name == "copy"
    assert copied.symmetry == symmetric_site.symmetry
    assert copied.grid is symmetric_site.grid
    assert copied.temperature == symmetric_site.temperature
    assert copied.bin_coords == symmetric_site.bin_coords

    for original_site, copied_site in zip(symmetric_site, copied):
        assert copied_site is not original_site
        assert copied_site.bin_coords == original_site.bin_coords
        assert copied_site.bin_coords is not original_site.bin_coords
        assert copied_site.site_counts_histogram is None
        assert copied_site.bulk_counts_histogram is None

    copied.get_site_list[0].bin_coords = {BinAddress(0, 0)}
    assert symmetric_site.get_site_list[0].bin_coords == {
        BinAddress(0, 7),
        BinAddress(1, 0),
    }


def test_rotate_bin_coords_rejects_site_number_outside_symmetry(symmetric_site):
    with pytest.raises(ValueError, match="greater than total symmetry"):
        symmetric_site._rotate_bin_coords({BinAddress(0, 0)}, 8, 4)
