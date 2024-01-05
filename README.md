# Polar_Binning_DeltaG
Descended from: https://github.com/BranniganLab/densitymap

Detailed protocol is in preparation.

# Basic usage:
1. Clone this repository
2. Install the DTA package with the following commands:
   ```
   cd Polar_Binning_DeltaG
   cd python
   #[Optional:] conda create -n DTA
   #[Optional:] conda activate DTA
   pip install .
   ```
3. Create or obtain unbiased coarse-grained (e.g. Martini) simulations of a membrane protein (only one protein currently possible without modifications)
4. Make sure the trajectory is whole (e.g. `gmx trjconv -f [inputname].xtc -s md.tpr -o [ouputname].xtc -ur compact -pbc mol -center`)
5. Load your trajectory into VMD
6. On first use for this protein: 

   a) write an assign_helices script based on one of the examples provided in ./TCL

   b) Test your helix assignment script from the tk console:
   
      ```> source <helix assignment script path>```

   c) In graphical representations, color the protein by occupancy and confirm that different colors correspond to different helices

7. Make a config file: 

    a) Select the closest config file from ./TCL/sample_configs and edit accordingly. 

    b) Test the selection strings that you used in the config file using the graphical representations window in VMD

8. From the tk console:

   ```> source TCL/polarDensity_for_DTA.tcl```

   ```> polarDensityBin <config file path>```  
9. Open the plot_enrichment.ipynb jupyter notebook using your method of choice (e.g. VSCode or a local host)
10. Update the paths and lipid names to point to your data and run
11. Select bins for the binding site (the show_bins notebook may be useful here)
12. Determine the accessible area for the site (many optional methods for this, get_accessible_area.ipynb describes one method)
13. Determine the expected mode of the bead probabilities in the bulk for a patch equal to the accessible area (e.g. using get_Npeak.ipynb)
14. Calculate the dG_bind using the get_probability_distributions.ipynb notebook
