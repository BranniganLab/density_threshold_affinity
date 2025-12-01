# -*- coding: utf-8 -*-
"""
Author:js2746

This is a pytest file for regression testing DTA.
"""

import pytest
from pathlib import Path

from DTA.Site import Site
from DTA.SymmetricSite import SymmetricSite
from DTA.SiteAcrossReplicas import SiteAcrossReplicas
from DTA.density import parse_tcl_dat_file, aggregate_density_enrichment_scores, load_replica_counts


def analyze_sample_outputs():
    # Make a SiteAcrossReplicas from sample outputs.
    root_path = Path("./test_files/sample_tcl_outputs/PC_CL").resolve()
    bulk_system_root_path = root_path.joinpath("../PC_CL_bulk").resolve()
    system_name = "CL"
    replicas = ["rep1", "rep2"]

    avg_enrichments = []
    for leaf in ["upp", "low"]:
        rep_list = []
        for rep in replicas:
            rep_path = root_path.joinpath(rep, f"{system_name}.{leaf}.avg.dat")
            rep_list.append(rep_path)
        avg_enrichment_over_reps, grid_dims = aggregate_density_enrichment_scores(rep_list)
        avg_enrichments.append(avg_enrichment_over_reps)

    site1 = Site(name="inner M1-M4", leaflet_id=2, temperature=320)
    site1.bin_coords = [(5, 13), (5, 14), (5, 15), (5, 16), (5, 17), (5, 18), (5, 19), (5, 20), (5, 21), (5, 22), (5, 23), (6, 18), (6, 19), (6, 20), (6, 21), (6, 22), (6, 23)]
    symm_site1 = SymmetricSite(symmetry=5, base_site=site1, Ntheta=grid_dims.Ntheta)
    replica_counts_list = load_replica_counts(root_path, replicas, system_name, site1.leaflet_id)
    symm_site_across_replicas_1 = SiteAcrossReplicas(replica_counts_list, base_site=symm_site1)

    bulk_counts_path = bulk_system_root_path.joinpath(f"{system_name}_counts_96.out")
    bulk_counts_list, _, _ = parse_tcl_dat_file(bulk_counts_path, bulk=True)
    symm_site_across_replicas_1.update_counts_histogram(bulk=True, counts_data=bulk_counts_list)

    return symm_site_across_replicas_1


def test_if_delta_G_matches():
    site = analyze_sample_outputs()
    assert round(site.dG, 2) == -1.45, "Delta G value no longer matches!"
    assert round(site.dG_std, 2) == 0.05, "Standard deviation no longer matches!"


@pytest.mark.xfail(strict=True)
def test_if_fail_appropriately():
    site = analyze_sample_outputs()
    assert site.dG == 0
