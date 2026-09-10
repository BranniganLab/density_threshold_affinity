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

    a) Copy one of the example configs from `./TCL/sample_configs/` and edit
    its values to match the loaded system and your desired parameters. The 
    examples are based on Martini ELIC and GPCR systems, but can be adapted.

    b) Test the atom selection strings that you use in the config file with
    VMD's graphical representations window.

7. From the tk console:

   ```> set dta_root /absolute/path/to/density_threshold_affinity```

   ```> source $dta_root/TCL/polarDensity_for_DTA.tcl```

   ```> polarDensityBin /absolute/path/to/your/config.tcl```
8. Open the Jupyter notebook that corresponds to your use case in density_threshold_affinity/python/ using your method of choice (e.g. VSCode or a local host)
9. Follow the prompts in the notebook

## Tcl configuration reference

The sample configuration files in `TCL/sample_configs/` provide runnable examples. The parameters they use are documented here so that the examples can remain concise.

### Leaflet assignment

- `leaflet_sorting_algorithm`: Selects the leaflet-assignment method.
  - `0`: Compare the specified head and tail selections.
  - `1`: Compare the default lipid termini.
  - `2`: Compare each lipid position with `leaflet_sorter_2_reference_sel`.
  - `3`: Compare each lipid with a local membrane midplane.
  - `4`: Assign every lipid to the upper leaflet.
- `leaflet_sorter_2_reference_sel`: VMD atom-selection string used as the reference for leaflet sorting algorithm 2. If set to `"none"`, the reference height is `z = 0`. Otherwise, the lipid center of mass is compared with the mass-weighted z center of this selection.
- `leaflet_reassign_interval`: Number of trajectory frames between recalculations of lipid leaflet identities. Between reassignment frames, the previous identities are reused.
- `restrict_leaflet_sorter_to_Rmax`: If enabled, restricts leaflet assignment to lipids within the configured radial analysis region.
- `headnames`: List of VMD atom-selection strings identifying the lipid head groups used by leaflet sorting algorithm 0. There must be one entry for each lipid selection in `atomsels`. If a selection contains multiple beads, their z coordinates are averaged before comparison with the corresponding tail selection.
- `tailnames`: List of VMD atom-selection strings identifying the lipid tail groups used by leaflet sorting algorithm 0. There must be one entry for each lipid selection in `atomsels`. If a selection contains multiple beads, their z coordinates are averaged before comparison with the corresponding head selection.

### Coordinate preparation

- `center_and_align`: If enabled, recenters and aligns the protein before density binning.
- `use_qwrap`: If enabled, uses the bundled qwrap library for coordinate wrapping. qwrap requires an orthorhombic unit cell. The bundled `qwrap.so` is loaded automatically.

### Helix assignment and protein geometry

- `helix_assignment_script`: Path to the system-specific Tcl script that assigns helix identities before analysis.
- `backbone_selstr`: VMD atom-selection string identifying the protein backbone beads used when calculating helix positions.
- `protein_selstr`: VMD atom-selection string identifying the protein atoms/beads used by the density-analysis workflow.
- `chainlist`: List of chain identifiers used when writing protein/inclusion coordinates. Each subunit or component that should be represented separately must have a distinct chain ID.
- `helixlist`: List of occupancy values assigned to the transmembrane helices by the helix-assignment script. These are occupancy values, not atom indices.
- `midplane_selstr`: VMD atom-selection string used to define the membrane midplane when protein/inclusion coordinates are divided between upper- and lower-leaflet outputs. The mass-weighted z center of the selection is calculated for each analyzed frame and then averaged over the trajectory. This parameter is separate from `leaflet_sorter_2_reference_sel`.

### Trajectory sampling

- `start_frame`: First trajectory frame included in the analysis.
- `end_frame`: Final trajectory frame included in the analysis. If omitted, the final loaded frame is used.
- `dt`: Interval, in frames, between trajectory frames included in the density calculation. `dt = 1` analyzes every frame. This is independent of `leaflet_reassign_interval`.

### Lipid selections and output

- `atomsels`: List of VMD atom-selection strings defining the lipid species or other membrane components for which densities should be calculated.
- `filename_stems`: List of filename stems used for outputs corresponding to `atomsels`. There must be one filename stem for each entry in `atomsels`.

### Polar grid

- `Rmin`: Minimum radial distance included in the polar density grid.
- `Rmax`: Maximum radial distance included in the polar density grid.
- `dr`: Radial bin width. `Rmax - Rmin` must be evenly divisible by `dr`.
- `Ntheta`: Number of angular bins in the polar density grid.
