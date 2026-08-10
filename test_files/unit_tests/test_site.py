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
    """Verify construction stores configuration but no results, because analysis must be explicit."""
    site = Site(name="site", grid=grid, leaflet_id=1, temperature=300)

    assert site.name == "site"
    assert site.grid is grid
    assert site.leaflet_id == 1
    assert site.temperature == 300
    assert site.bin_coords is None
    assert site.site_counts_histogram is None
    assert site.bulk_counts_histogram is None


def test_init_rejects_non_grid():
    """Verify a non-PolarBinGrid is rejected, because Site geometry depends on the grid API."""
    with pytest.raises(TypeError, match="grid must be a PolarBinGrid"):
        Site("site", object(), 1, 300)


@pytest.mark.parametrize("leaflet_id", [0, 3, "upper", None])
def test_init_rejects_invalid_leaflet_id(grid, leaflet_id):
    """Verify unsupported leaflet identifiers fail, preventing ambiguous leaflet assignments."""
    with pytest.raises(ValueError, match="leaflet_id must be 1 or 2"):
        Site("site", grid, leaflet_id, 300)


@pytest.mark.parametrize("temperature", [0, -1, -273.15])
def test_init_rejects_nonpositive_temperature(grid, temperature):
    """Verify nonpositive Kelvin temperatures fail because they are physically invalid."""
    with pytest.raises(ValueError, match="temperature must be positive"):
        Site("site", grid, 1, temperature)


@pytest.mark.parametrize("container_type", [list, tuple, set])
def test_bin_coords_accepts_supported_containers_and_converts_tuples(grid, container_type):
    """Verify supported collections normalize to BinAddress sets for consistent downstream use."""
    site = Site("site", grid, 1, 300)
    site.bin_coords = container_type([(0, 1), BinAddress(1, 2), (0, 1)])

    assert site.bin_coords == {BinAddress(0, 1), BinAddress(1, 2)}
    assert all(isinstance(address, BinAddress) for address in site.bin_coords)


def test_bin_coords_rejects_invalid_container(site):
    """Verify unsupported containers fail clearly instead of being misinterpreted as bin addresses."""
    with pytest.raises(TypeError, match="list, tuple, set, or frozenset"):
        site.bin_coords = np.array([[0, 0]])


@pytest.mark.parametrize("empty_container", [[], (), set()])
def test_bin_coords_rejects_empty_collections(grid, empty_container):
    """Verify a Site cannot be defined without bins because it could not measure occupancy."""
    site = Site("site", grid, 1, 300)

    with pytest.raises(ValueError, match="bin_addresses cannot be empty"):
        site.bin_coords = empty_container


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
    """Verify invalid radial and angular indices fail before they can index count arrays."""
    site = Site("site", grid, 1, 300)

    with pytest.raises(IndexError, match=message):
        site.bin_coords = [address]


def test_update_bulk_counts_histogram(site):
    """Verify raw bulk counts are histogrammed and expose the modal occupancy as n_peak."""
    site.update_bulk_counts_histogram(np.array([0, 1, 1, 3]))

    np.testing.assert_array_equal(site.bulk_counts_histogram, np.array([1, 2, 0, 1]))
    assert site.n_peak == 1


def test_n_peak_requires_bulk_histogram(site):
    """Verify n_peak fails before bulk data are supplied, because its value is bulk-derived."""
    with pytest.raises(AttributeError, match="update the bulk counts histogram"):
        _ = site.n_peak


@pytest.mark.parametrize(
    "update_method_name",
    ["update_site_counts_histogram", "update_bulk_counts_histogram"],
)
def test_histogram_updates_require_defined_bin_coords(grid, update_method_name):
    """Verify histogram updates fail while the Site geometry remains undefined."""
    site = Site("site", grid, 1, 300)
    counts = np.zeros((3, 2, 8), dtype=int)
    if update_method_name == "update_bulk_counts_histogram":
        counts = np.zeros(3, dtype=int)

    with pytest.raises(RuntimeError, match="before bin_coords have been defined"):
        getattr(site, update_method_name)(counts)


def test_update_site_counts_histogram_sums_selected_bins_over_time(site):
    """Verify selected-bin counts are summed per frame and histogrammed for affinity analysis."""
    counts = np.zeros((4, 2, 8), dtype=float)
    counts[:, 0, 0] = [0.0, 1.0, 1.0, 2.0]
    counts[:, 1, 1] = [1.0, 0.0, 2.0, 1.0]

    site.update_site_counts_histogram(counts)

    # The input is converted to integers before the selected bins are summed.
    np.testing.assert_array_equal(site.site_counts_histogram, np.array([0, 2, 0, 2]))


@pytest.mark.parametrize(
    ("update_method_name", "counts"),
    [
        ("update_site_counts_histogram", [[[0]]]),
        ("update_bulk_counts_histogram", [0, 1]),
    ],
)
def test_histogram_updates_reject_non_arrays(site, update_method_name, counts):
    """Verify both update paths require ndarrays so dtype and dimensions are unambiguous."""
    with pytest.raises(TypeError, match="must be provided as a NumPy ndarray"):
        getattr(site, update_method_name)(counts)


def test_update_bulk_counts_histogram_rejects_non_1d_array(site):
    """Verify bulk input must be one-dimensional because it represents one count per frame."""
    with pytest.raises(ValueError, match="counts_data must be 1D"):
        site.update_bulk_counts_histogram(np.zeros((2, 2), dtype=int))


def test_update_site_counts_histogram_rejects_non_3d_array(site):
    """Verify site input must have time, radial, and angular axes to locate selected bins."""
    with pytest.raises(ValueError, match="counts_data must be 3D"):
        site.update_site_counts_histogram(np.zeros((2, 8), dtype=int))


def test_update_site_counts_histogram_rejects_wrong_grid_shape(site):
    """Verify count-grid dimensions match the Site grid so bin addresses retain their meaning."""
    with pytest.raises(ValueError, match="counts_data is the wrong shape"):
        site.update_site_counts_histogram(np.zeros((3, 3, 8), dtype=int))


@pytest.mark.parametrize(
    ("update_method_name", "counts"),
    [
        ("update_site_counts_histogram", np.full((2, 2, 8), "one")),
        ("update_bulk_counts_histogram", np.array(["zero", "one"])),
        ("update_bulk_counts_histogram", np.array([True, False])),
        ("update_bulk_counts_histogram", np.array([0 + 0j, 1 + 0j])),
    ],
)
def test_histogram_updates_reject_nonnumeric_arrays(site, update_method_name, counts):
    """Verify text, Boolean, and complex arrays fail because counts must be real numbers."""
    with pytest.raises(TypeError, match="must contain real numeric values"):
        getattr(site, update_method_name)(counts)


@pytest.mark.parametrize("bad_value", [np.nan, np.inf, -np.inf])
@pytest.mark.parametrize(
    ("update_method_name", "shape"),
    [
        ("update_site_counts_histogram", (2, 2, 8)),
        ("update_bulk_counts_histogram", (2,)),
    ],
)
def test_histogram_updates_reject_nonfinite_counts(
    site,
    update_method_name,
    shape,
    bad_value,
):
    """Verify NaN and infinite values fail because they cannot represent molecule counts."""
    counts = np.zeros(shape, dtype=float)
    counts.flat[0] = bad_value

    with pytest.raises(ValueError, match="must contain only finite values"):
        getattr(site, update_method_name)(counts)


@pytest.mark.parametrize(
    ("update_method_name", "counts"),
    [
        ("update_site_counts_histogram", np.full((2, 2, 8), 1.5)),
        ("update_bulk_counts_histogram", np.array([0.0, 1.5])),
    ],
)
def test_histogram_updates_reject_fractional_counts(site, update_method_name, counts):
    """Verify fractional values fail rather than being silently truncated to integers."""
    with pytest.raises(ValueError, match="must contain integer-valued counts"):
        getattr(site, update_method_name)(counts)


@pytest.mark.parametrize(
    ("update_method_name", "counts"),
    [
        ("update_site_counts_histogram", -np.ones((2, 2, 8), dtype=int)),
        ("update_bulk_counts_histogram", np.array([0, -1], dtype=int)),
    ],
)
def test_histogram_updates_reject_negative_counts(site, update_method_name, counts):
    """Verify negative values fail because molecule counts cannot be below zero."""
    with pytest.raises(ValueError, match="cannot contain negative counts"):
        getattr(site, update_method_name)(counts)


def test_histogram_updates_accept_integer_valued_float_arrays(site):
    """Verify integral floats are accepted and normalized when parsers return float arrays."""
    site_counts = np.zeros((3, 2, 8), dtype=float)
    site_counts[:, 0, 0] = [0.0, 1.0, 2.0]
    bulk_counts = np.array([0.0, 1.0, 1.0])

    site.update_site_counts_histogram(site_counts)
    site.update_bulk_counts_histogram(bulk_counts)

    np.testing.assert_array_equal(site.site_counts_histogram, np.array([1, 1, 1]))
    np.testing.assert_array_equal(site.bulk_counts_histogram, np.array([1, 2]))


def test_dg_requires_site_histogram(site):
    """Verify dG cannot be evaluated without site occupancy statistics."""
    site.update_bulk_counts_histogram(np.array([0, 1, 1, 2]))

    with pytest.raises(AssertionError, match="site counts histogram"):
        _ = site.dG


def test_dg_requires_bulk_histogram(site):
    """Verify dG cannot be evaluated without the bulk reference correction."""
    site.update_site_counts_histogram(np.zeros((3, 2, 8), dtype=int))

    with pytest.raises(AssertionError, match=r"add bulk\s+counts"):
        _ = site.dG


def test_dg_subtracts_bulk_reference_from_site_free_energy(site):
    """Verify dG is the site free energy minus the bulk-reference free energy."""
    counts = np.zeros((5, 2, 8), dtype=int)
    counts[:, 0, 0] = [0, 1, 2, 2, 3]
    site.update_site_counts_histogram(counts)
    site.update_bulk_counts_histogram(np.array([0, 1, 1, 2, 2]))

    expected = calculate_dG(site.site_counts_histogram, 1, 320)
    expected -= calculate_dG(site.bulk_counts_histogram, 1, 320)
    assert site.dG == pytest.approx(expected)


def test_calculate_geometric_area_sums_selected_bin_areas(site):
    """Verify geometric area sums every selected polar bin rather than using bin count alone."""
    expected = site.grid.calc_bin_area(BinAddress(0, 0))
    expected += site.grid.calc_bin_area(BinAddress(1, 1))

    assert site.calculate_geometric_area() == pytest.approx(expected)


def test_predict_accessible_area_uses_histogram_modes(site):
    """Verify the default area estimate scales bulk area by modal site-to-bulk occupancy."""
    counts = np.zeros((5, 2, 8), dtype=int)
    counts[:, 0, 0] = [1, 2, 2, 2, 3]
    site.update_site_counts_histogram(counts)
    site.update_bulk_counts_histogram(np.array([1, 1, 2, 2, 2]))

    assert site.predict_accessible_area(100.0) == pytest.approx(100.0)


def test_predict_accessible_area_can_use_histogram_means(site):
    """Verify callers can use mean occupancies when modes discard useful distribution detail."""
    counts = np.zeros((5, 2, 8), dtype=int)
    counts[:, 0, 0] = [1, 2, 2, 2, 3]
    site.update_site_counts_histogram(counts)
    site.update_bulk_counts_histogram(np.array([1, 1, 2, 2, 2]))

    # Site mean = 2.0 and bulk mean = 1.6.
    assert site.predict_accessible_area(100.0, mode=False) == pytest.approx(125.0)


def test_copy_preserves_definition_but_clears_histograms(site):
    """Verify copying reuses a Site definition without carrying ligand-specific analysis results."""
    site.update_site_counts_histogram(np.zeros((3, 2, 8), dtype=int))
    site.update_bulk_counts_histogram(np.array([0, 1, 1]))

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
    """Verify a partially defined Site can be copied without inventing a bin selection."""
    site = Site("site", grid, 1, 300)

    copied = site.copy("copy")

    assert copied.bin_coords is None
