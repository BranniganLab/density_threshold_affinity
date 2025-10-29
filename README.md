# Density Threshold Affinity

Detailed protocol is published! 

Check it out here: https://doi.org/10.1016/bs.mie.2024.03.008 


# Basic usage:
1. Clone this repository
2. Install the DTA package with the following commands:
   ```
   cd density_threshold_affinity
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

   ```> source ./TCL/polarDensity_for_DTA.tcl```

   ```> polarDensity <config file path>```  
9. Open the jupyter notebook that corresponds to your use case using your method of choice (e.g. VSCode or a local host)
10. Follow the prompts in the notebook
