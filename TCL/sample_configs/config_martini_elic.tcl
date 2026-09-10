# Configuration example based on a Martini ELIC system.
#
# Edit the atom selections, paths, chain names, and helix assignments to match
# the loaded molecule before running polarDensityBin. See the Tcl configuration
# reference in README.md for parameter documentation.

# Leaflet assignment
set leaflet_sorting_algorithm 1
set leaflet_sorter_2_reference_sel "none"
set leaflet_reassign_interval 1
set restrict_leaflet_sorter_to_Rmax 0

# Coordinate preparation
set center_and_align 0
set use_qwrap 0

# Helix assignment script
set helix_assignment_script "assign_helices_ELIC_general.tcl"

# Frames and trajectory sampling
set start_frame 0
set dt 1

# Protein selections and helix assignment
set backbone_selstr "name BB"
set protein_selstr "name BB SC1 to SC4"
# The bundled ELIC helix-assignment example assigns the five subunits chains A-E.
set chainlist [list A B C D E]
set helixlist [list 1 2 3 4]
set midplane_selstr "occupancy 1 to 4"

# Lipid selections and output names
set atomsels [list "resname POPG"]
set filename_stems [list "POPG"]
set headnames [list "name PO4"]
set tailnames [list "name C4"]

# Polar grid
set Rmax 20.
set Rmin 0.
set dr 5.
set Ntheta 50
