#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for the Site class."""

import numpy as np
import pytest

from dta.Site import Site
from dta.bin_logic import BinAddress, PolarBinGrid
from dta.utils import calculate_dG


@pytest.fixture
def grid():
    """Return a small polar grid for Site tests."""
    return PolarBinGrid(r_min=0.0, r_max=2.0, n_r=2, n_theta=8)


@pytest.fixture
def site(grid):
    """Return a Site containing two bins."""
    result = Site(name="site", grid=grid, leaflet_id=2, temperature=320)
    result.bin_coords = {BinAddress(0, 0), BinAddress(1, 1)}
    return result


def test_init_sets_attributes_and_empty_analysis_state(grid):
    site = Site(name="site", grid=grid, leaflet_id=1, temperature=300)

    assert site.name == "site"
    assert site.grid is grid
    assert site.leaflet_id == 1
    assert site.temperature == 300
    assert site.bin_coords is None
    assert site.site_counts_histogram is None
    assert site.bulk_counts_histogram is None


def test_init_rejects_non_grid():
    with pytest.raises(TypeError, match="grid must be a PolarBinGrid"):
        Site("site", object(), 1, 300)


@pytest.mark.parametrize("leaflet_id", [0, 3, "upper", None])
def test_init_rejects_invalid_leaflet_id(grid, leaflet_id):
    with pytest.raises(ValueError, match="leaflet_id must be 1 or 2"):
        Site("site", grid, leaflet_id, 300)


@pytest.mark.parametrize("container_type", [list, tuple, set])
def test_bin_coords_accepts_supported_containers_and_converts_tuples(grid, container_type):
    site = Site("site", grid, 1, 300)
    site.bin_coords = container_type([(0, 1), BinAddress(1, 2), (0, 1)])

    assert site.bin_coords == {BinAddress(0, 1), BinAddress(1, 2)}
    assert all(isinstance(address, BinAddress) for address in site.bin_coords)


def test_bin_coords_rejects_invalid_container(site):
    with pytest.raises(TypeError, match="list, tuple, or set"):
        site.bin_coords = np.array([[0, 0]])


@pytest.mark.parametrize(
    ("address", "message"),
    [
        ((-1, 0), "Radial bin -1 out of range"),
        ((2, 0), "Radial bin 2 out of range"),
        ((0, -1), "Angular bin -1 out of range"),
        ((0, 8), "Angular bin 8 out of range"),
    ],
)
def test_bin_coords_rejects_out_of_range_addresses(grid, address, message):
    site = Site("site", grid, 1, 300)

    with pytest.raises(IndexError, match=message):
        site.bin_coords = [address]


def test_update_bulk_counts_histogram(site):
    site.update_counts_histogram(bulk=True, counts_data=np.array([0, 1, 1, 3]))

    np.testing.assert_array_equal(site.bulk_counts_histogram, np.array([1, 2, 0, 1]))
    assert site.n_peak == 1


def test_n_peak_requires_bulk_histogram(site):
    with pytest.raises(AttributeError, match="update the bulk counts histogram"):
        _ = site.n_peak


def test_update_site_counts_histogram_sums_selected_bins_over_time(site):
    counts = np.zeros((4, 2, 8), dtype=float)
    counts[:, 0, 0] = [0.0, 1.0, 1.9, 2.0]
    counts[:, 1, 1] = [1.0, 0.0, 2.1, 1.0]

    site.update_counts_histogram(bulk=False, counts_data=counts)

    # The input is converted to integers before the selected bins are summed.
    np.testing.assert_array_equal(site.site_counts_histogram, np.array([0, 2, 0, 2]))


def test_update_counts_histogram_rejects_non_array(site):
    with pytest.raises(TypeError, match="ndarray not supplied"):
        site.update_counts_histogram(bulk=True, counts_data=[0, 1])


def test_update_bulk_counts_histogram_rejects_non_1d_array(site):
    with pytest.raises(ValueError, match="Bulk counts data is not in the right format"):
        site.update_counts_histogram(bulk=True, counts_data=np.zeros((2, 2), dtype=int))


def test_update_site_counts_histogram_rejects_non_3d_array(site):
    with pytest.raises(ValueError, match="Counts data is not in the right format"):
        site.update_counts_histogram(bulk=False, counts_data=np.zeros((2, 8), dtype=int))


def test_update_site_counts_histogram_rejects_wrong_grid_shape(site):
    with pytest.raises(ValueError, match="counts_data is the wrong shape"):
        site.update_counts_histogram(bulk=False, counts_data=np.zeros((3, 3, 8), dtype=int))


def test_dg_requires_site_histogram(site):
    site.update_counts_histogram(bulk=True, counts_data=np.array([0, 1, 1, 2]))

    with pytest.raises(AssertionError, match="site counts histogram"):
        _ = site.dG


def test_dg_requires_bulk_histogram(site):
    site.update_counts_histogram(bulk=False, counts_data=np.zeros((3, 2, 8), dtype=int))

    with pytest.raises(AssertionError, match=r"add bulk\s+counts"):
        _ = site.dG


def test_dg_subtracts_bulk_reference_from_site_free_energy(site):
    counts = np.zeros((5, 2, 8), dtype=int)
    counts[:, 0, 0] = [0, 1, 2, 2, 3]
    site.update_counts_histogram(bulk=False, counts_data=counts)
    site.update_counts_histogram(bulk=True, counts_data=np.array([0, 1, 1, 2, 2]))

    expected = calculate_dG(site.site_counts_histogram, 1, 320)
    expected -= calculate_dG(site.bulk_counts_histogram, 1, 320)
    assert site.dG == pytest.approx(expected)


def test_calculate_geometric_area_sums_selected_bin_areas(site):
    expected = site.grid.calc_bin_area(BinAddress(0, 0))
    expected += site.grid.calc_bin_area(BinAddress(1, 1))

    assert site.calculate_geometric_area() == pytest.approx(expected)


def test_predict_accessible_area_uses_histogram_modes(site):
    counts = np.zeros((5, 2, 8), dtype=int)
    counts[:, 0, 0] = [1, 2, 2, 2, 3]
    site.update_counts_histogram(bulk=False, counts_data=counts)
    site.update_counts_histogram(bulk=True, counts_data=np.array([1, 1, 2, 2, 2]))

    assert site.predict_accessible_area(100.0) == pytest.approx(100.0)


def test_predict_accessible_area_can_use_histogram_means(site):
    counts = np.zeros((5, 2, 8), dtype=int)
    counts[:, 0, 0] = [1, 2, 2, 2, 3]
    site.update_counts_histogram(bulk=False, counts_data=counts)
    site.update_counts_histogram(bulk=True, counts_data=np.array([1, 1, 2, 2, 2]))

    # Site mean = 2.0 and bulk mean = 1.6.
    assert site.predict_accessible_area(100.0, mode=False) == pytest.approx(125.0)


def test_copy_preserves_definition_but_clears_histograms(site):
    site.update_counts_histogram(bulk=False, counts_data=np.zeros((3, 2, 8), dtype=int))
    site.update_counts_histogram(bulk=True, counts_data=np.array([0, 1, 1]))

    copied = site.copy("copy")

    assert copied is not site
    assert copied.name == "copy"
    assert copied.grid is site.grid
    assert copied.leaflet_id == site.leaflet_id
    assert copied.temperature == site.temperature
    assert copied.bin_coords == site.bin_coords
    assert copied.bin_coords is not site.bin_coords
    assert copied.site_counts_histogram is None
    assert copied.bulk_counts_histogram is None

    copied.bin_coords = {BinAddress(0, 2)}
    assert site.bin_coords == {BinAddress(0, 0), BinAddress(1, 1)}


def test_copy_preserves_undefined_bin_coords(grid):
    site = Site("site", grid, 1, 300)

    copied = site.copy("copy")

    assert copied.bin_coords is None
