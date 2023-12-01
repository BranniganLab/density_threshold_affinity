# Polar_Binning_DeltaG
Descended from: https://github.com/BranniganLab/densitymap

Detailed protocol is in preparation.

# Basic usage:
1. Clone this repository
2. Create or obtain unbiased coarse-grained (e.g. Martini) simulations of a membrane protein (only one protein currently possible without modifications)
3. Make sure the trajectory is whole (e.g. `gmx trjconv -f [inputname].xtc -s md.tpr -o [ouputname].xtc -ur compact -pbc mol -center`)
4. Write an assign_helices script based on one of the examples provided in ./TCL
5. Update the do_polar_density_for_DTA.tcl script in ./TCL according to the documentation therein
6. Load your trajectory into VMD
7. Source do_polar_density_for_DTC.tcl from the tkconsole
8. Open the plot_enrichment.ipynb jupyter notebook using your method of choice (e.g. VSCode or a local host)
9. Update the paths and lipid names to point to your data and run
10. Select bins for the binding site (the show_bins notebook may be useful here)
11. Determine the accessible area for the site (many optional methods for this, get_accessible_area.ipynb describes one method)
12. Determine the expected mode of the bead probabilities in the bulk for a patch equal to the accessible area (e.g. using get_Npeak.ipynb)
13. Calculate the dG_bind using the get_probability_distributions.ipynb notebook
