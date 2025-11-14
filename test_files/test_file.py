# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pytest
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

from DTA.utils import load_inclusion_helices
from DTA.density import parse_tcl_dat_file, aggregate_density_enrichment_scores, load_replica_counts, valid_Dimensions
from DTA.plotting import make_density_enrichment_heatmap, make_custom_colormap, plot_histogram, outline_site, plot_titration_curve
from DTA.Site import Site
from DTA.SymmetricSite import SymmetricSite
from DTA.SiteAcrossReplicas import SiteAcrossReplicas

colormap = make_custom_colormap()
root_path = Path("./test_files/sample_tcl_outputs/PC_CL").resolve() # The directory containing your replica subdirectories
DPPC_root_path = root_path.joinpath("../DPPC").resolve() # The directory containing the outputs of PolarDensityBin from your DPPC + protein simulation
DPPC_bulk_root_path = root_path.joinpath("../DPPC_bulk/").resolve() # The directory that will contain the outputs of do_get_counts.tcl from your DPPC bulk simulation
bulk_system_root_path = root_path.joinpath("../PC_CL_bulk").resolve() # The directory that will contain the outputs of do_get_counts.tcl from your bulk simulation
system_name = "CL" # The filestem that PolarDensityBin used.
replicas = ["rep1", "rep2"] # names of replica subdirectories, located in the "root_path" directory
helix_definitions = root_path.joinpath(replicas[0]) #where are the coordinates for the transmembrane helices?
max_enrichment = 5 # how high do you want your heat map to go?

print(str(root_path))

helices = load_inclusion_helices(helix_definitions)
avg_enrichments = []
for leaf in ["upp", "low"]:
    rep_list = []
    for rep in replicas:
        rep_path = root_path.joinpath(rep, f"{system_name}.{leaf}.avg.dat")
        rep_list.append(rep_path)
    avg_enrichment_over_reps, grid_dims = aggregate_density_enrichment_scores(rep_list)
    avg_enrichments.append(avg_enrichment_over_reps)

print(helices)
print(avg_enrichments)