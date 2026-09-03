# Density Threshold Affinity

Detailed protocol is published! 

Check it out here: https://doi.org/10.1016/bs.mie.2024.03.008 


# Basic usage:
1. Clone this repository
2. Install the DTA package with the following commands:
   ```
   cd density_threshold_affinity
   cd python
   pip install .
   ```
3. Create or obtain coarse-grained (e.g. Martini) simulations of a membrane protein (only one protein currently possible without modifications)
4. Load your trajectory into VMD
5. On first use for this protein: 

   a) write an assign_helices script based on one of the examples provided in ./TCL/sample_helix_assignment_scripts/

   b) Test your helix assignment script from the tk console:
   
      ```> source <helix assignment script path>```

   c) In graphical representations, color the protein by occupancy and confirm that different colors correspond to different helices

6. Make a config file:

    a) Select a sample config from `./TCL/sample_configs/` based on the
    system being analyzed:

    - Use `config_martini_elic.tcl` for an ELIC-like pentameric system.
    - Use `config_martini_gpcr.tcl` for a GPCR-like system.

    Edit the selected config to match the loaded system. In particular, check
    the protein selections, chain and helix lists, helix-assignment script,
    lipid selections, and utility paths. The sample configs document every
    supported parameter and provide working examples for these two system
    types.

    b) Test the atom selection strings that you use in the config file with
    VMD's graphical representations window.

   The configuration uses these groups of parameters:

   - `leaflet_sorting_algorithm`, `leaflet_sorter_2_reference_sel`,
     `leaflet_reassign_interval`, and `restrict_leaflet_sorter_to_Rmax`
     control leaflet assignment.
   - `center_and_align`, `use_qwrap`, and `utils` control coordinate
     preparation and the VMD utility scripts.
   - `start_frame`, `end_frame`, and `dt` select the trajectory frames.
   - `backbone_selstr`, `protein_selstr`, `chainlist`, `helixlist`,
     `helix_assignment_script`, and `midplane_selstr` describe the protein.
     The helix-assignment script assigns occupancy values to the transmembrane
     helices; `helixlist` is the list of those occupancy values, not atom
     indices. `chainlist` contains the chain identifiers used when writing
     protein coordinates, so a system with one chain containing all subunits
     should use a one-item chain list.
   - `atomsels` contains one lipid atom selection per species and
     `filename_stems` contains the matching output name for each selection.
     `headnames` and `tailnames` are used by leaflet sorting algorithm 0.
     Each is an atom-selection string; multiple beads can be selected, for
     example `name C4 C5 C6`. Their centers of mass are compared to assign
     the leaflet.
   - `Rmin`, `Rmax`, `dr`, and `Ntheta` define the polar grid. `Rmax - Rmin`
     must be evenly divisible by `dr`.

7. From the tk console:

   ```> source ./TCL/polarDensity_for_DTA.tcl```

   ```> polarDensityBin <config file path>```
8. Open the Jupyter notebook that corresponds to your use case in density_threshold_affinity/python/ using your method of choice (e.g. VSCode or a local host)
9. Follow the prompts in the notebook
