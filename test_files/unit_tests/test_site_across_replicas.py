#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for the SiteAcrossReplicas class."""

import numpy as np
import pytest

from dta.Site import Site
from dta.SiteAcrossReplicas import SiteAcrossReplicas
from dta.SymmetricSite import SymmetricSite
from dta.bin_logic import BinAddress, PolarBinGrid
from dta.utils import calculate_dG


@pytest.fixture
def grid():
    """Return a small grid that also supports fourfold symmetric sites."""
    return PolarBinGrid(r_min=0.0, r_max=2.0, n_r=2, n_theta=8)


@pytest.fixture
def base_site(grid):
    """Return a fully defined Site that can be copied across replicas."""
    site = Site(name="binding site", grid=grid, leaflet_id=2, temperature=320)
    site.bin_coords = {BinAddress(0, 0), BinAddress(1, 1)}
    return site


@pytest.fixture
def replica_counts():
    """Return two replicas with distinct but finite site-occupancy distributions."""
    replica_1 = np.zeros((5, 2, 8), dtype=int)
    replica_1[:, 0, 0] = [0, 1, 1, 2, 2]

    replica_2 = np.zeros((5, 2, 8), dtype=int)
    replica_2[:, 0, 0] = [0, 1, 2, 2, 3]

    return [replica_1, replica_2]


@pytest.fixture
def site_across_replicas(base_site, replica_counts):
    """Return a SiteAcrossReplicas populated with site counts but no bulk counts."""
    return SiteAcrossReplicas(replica_list=replica_counts, base_site=base_site)


@pytest.fixture
def populated_site_across_replicas(site_across_replicas):
    """Return a SiteAcrossReplicas with the shared bulk histogram populated."""
    site_across_replicas.update_bulk_counts_histogram(
        np.array([0, 1, 1, 2, 2]),
    )
    return site_across_replicas


def test_init_creates_one_named_site_per_replica_and_inherits_configuration(
    base_site,
    replica_counts,
):
    """Verify construction preserves the site definition while naming one copy per replica."""
    across_replicas = SiteAcrossReplicas(replica_counts, base_site)

    assert across_replicas.name == "binding site"
    assert across_replicas.grid is base_site.grid
    assert across_replicas.temperature == base_site.temperature
    assert [site.name for site in across_replicas] == [
        "binding site_rep1",
        "binding site_rep2",
    ]
    assert all(site is not base_site for site in across_replicas)
    assert all(site.bin_coords == base_site.bin_coords for site in across_replicas)
    assert base_site.site_counts_histogram is None


def test_init_loads_each_replica_into_its_corresponding_site(
    site_across_replicas,
):
    """Verify replica-specific count arrays produce distinct constituent histograms."""
    first, second = site_across_replicas.get_site_list

    np.testing.assert_array_equal(first.site_counts_histogram, np.array([1, 2, 2]))
    np.testing.assert_array_equal(second.site_counts_histogram, np.array([1, 1, 2, 1]))


def test_iter_yields_constituent_sites_in_replica_order(site_across_replicas):
    """Verify iteration retains replica order so results remain associated with input replicas."""
    assert list(site_across_replicas) == site_across_replicas.get_site_list
    assert all(isinstance(site, Site) for site in site_across_replicas)


def test_init_supports_symmetric_site_definitions(grid, replica_counts):
    """Verify a SymmetricSite can be copied per replica without losing its rotations."""
    base_site = Site("symmetric binding site", grid, leaflet_id=2, temperature=320)
    base_site.bin_coords = {BinAddress(0, 0)}
    symmetric_site = SymmetricSite(symmetry=4, base_site=base_site)

    across_replicas = SiteAcrossReplicas(replica_counts, symmetric_site)

    assert across_replicas.name == "symmetric binding site"
    assert len(across_replicas.get_site_list) == 2
    assert all(isinstance(site, SymmetricSite) for site in across_replicas)
    assert all(site.symmetry == 4 for site in across_replicas)
    assert [site.name for site in across_replicas] == [
        "symmetric binding site_rep1",
        "symmetric binding site_rep2",
    ]
    assert all(site.bin_coords == symmetric_site.bin_coords for site in across_replicas)


def test_init_rejects_site_without_bin_coordinates(grid, replica_counts):
    """Verify an undefined Site fails because its replica counts cannot be spatially selected."""
    base_site = Site("undefined", grid, leaflet_id=1, temperature=300)

    with pytest.raises(ValueError, match="base_site needs to be fully defined"):
        SiteAcrossReplicas(replica_counts, base_site)


@pytest.mark.parametrize("base_site", [None, object(), "site"])
def test_init_rejects_unsupported_base_site_types(replica_counts, base_site):
    """Verify only Site and SymmetricSite definitions are accepted for replica construction."""
    with pytest.raises(ValueError, match="base_site must be a Site or SymmetricSite"):
        SiteAcrossReplicas(replica_counts, base_site)


def test_init_rejects_non_list_replica_collection(base_site, replica_counts):
    """Verify replica arrays must be supplied in a list so their ordering is explicit."""
    with pytest.raises(TypeError, match="replica_list must be a list"):
        SiteAcrossReplicas(tuple(replica_counts), base_site)


def test_init_rejects_empty_replica_list(base_site):
    """Verify construction requires at least one replica so aggregation remains meaningful."""
    with pytest.raises(ValueError, match="replica_list cannot be empty"):
        SiteAcrossReplicas([], base_site)


def test_init_propagates_invalid_replica_array_shape(base_site):
    """Verify malformed replica arrays fail rather than creating unusable constituent Sites."""
    with pytest.raises(ValueError, match="counts_data must be 3D"):
        SiteAcrossReplicas([np.zeros((2, 8), dtype=int)], base_site)


@pytest.mark.parametrize(
    ("replica_counts", "error_type", "message"),
    [
        ([[[0]]], TypeError, "must be provided as a NumPy ndarray"),
        ([np.full((2, 2, 8), "one")], TypeError, "must contain real numeric values"),
        ([np.full((2, 2, 8), np.nan)], ValueError, "must contain only finite values"),
        ([np.full((2, 2, 8), 1.5)], ValueError, "must contain integer-valued counts"),
        ([-np.ones((2, 2, 8), dtype=int)], ValueError, "cannot contain negative counts"),
    ],
)
def test_init_propagates_replica_count_validation(
    base_site,
    replica_counts,
    error_type,
    message,
):
    """Verify each replica is subject to Site count validation during construction."""
    with pytest.raises(error_type, match=message):
        SiteAcrossReplicas(replica_counts, base_site)


def test_site_counts_histogram_aggregates_all_replicas(site_across_replicas):
    """Verify the reported histogram sums replica histograms, including unequal lengths."""
    np.testing.assert_array_equal(
        site_across_replicas.site_counts_histogram,
        np.array([2.0, 3.0, 4.0, 1.0]),
    )


def test_update_bulk_counts_histogram_updates_every_replica(site_across_replicas):
    """Verify a shared bulk dataset is propagated to every replica for consistent correction."""
    bulk_counts = np.array([0, 1, 1, 2, 2])

    site_across_replicas.update_bulk_counts_histogram(bulk_counts)

    for site in site_across_replicas:
        np.testing.assert_array_equal(site.bulk_counts_histogram, np.array([1, 2, 2]))


def test_update_site_counts_histogram_replaces_each_replica_with_shared_data(
    site_across_replicas,
):
    """Verify a non-bulk update is deliberately applied to every constituent replica Site."""
    replacement = np.zeros((4, 2, 8), dtype=int)
    replacement[:, 0, 0] = [0, 1, 1, 2]

    site_across_replicas.update_site_counts_histogram(replacement)

    for site in site_across_replicas:
        np.testing.assert_array_equal(site.site_counts_histogram, np.array([1, 2, 1]))


def test_bulk_counts_histogram_and_n_peak_use_shared_bulk_distribution(
    populated_site_across_replicas,
):
    """Verify the common bulk histogram is exposed and supplies the occupancy threshold."""
    np.testing.assert_array_equal(
        populated_site_across_replicas.bulk_counts_histogram,
        np.array([1, 2, 2]),
    )
    assert populated_site_across_replicas.n_peak == 1


def test_dg_uses_aggregated_replica_and_bulk_histograms(
    populated_site_across_replicas,
):
    """Verify overall affinity compares pooled replica occupancy with the bulk reference."""
    expected = calculate_dG(
        populated_site_across_replicas.site_counts_histogram,
        n_peak=1,
        temperature=320,
    )
    expected -= calculate_dG(
        populated_site_across_replicas.bulk_counts_histogram,
        n_peak=1,
        temperature=320,
    )

    assert populated_site_across_replicas.dG == pytest.approx(expected)


def test_dg_std_is_standard_deviation_across_replica_affinities(
    populated_site_across_replicas,
):
    """Verify dG_std reports between-replica affinity variation as a convergence diagnostic."""
    expected = np.std(
        np.array([site.dG for site in populated_site_across_replicas]),
    )

    assert populated_site_across_replicas.dG_std == pytest.approx(expected)


def test_predict_accessible_area_uses_histogram_modes(
    populated_site_across_replicas,
):
    """Verify default area prediction scales bulk area by pooled-site and bulk modes."""
    # The pooled site mode is 2, while the bulk mode is 1.
    assert populated_site_across_replicas.predict_accessible_area(100.0) == pytest.approx(200.0)


def test_predict_accessible_area_can_use_histogram_means(
    populated_site_across_replicas,
):
    """Verify mean-based prediction uses the full pooled replica occupancy distribution."""
    # The pooled site mean is 1.4, while the bulk mean is 1.2.
    expected = 100.0 * 1.4 / 1.2
    assert populated_site_across_replicas.predict_accessible_area(
        100.0,
        mode=False,
    ) == pytest.approx(expected)
