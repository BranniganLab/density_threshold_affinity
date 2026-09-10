# Configuration example based on a Martini GPCR system.
#
# Edit the atom selections, paths, chain names, and helix assignments to match
# the loaded molecule before running polarDensityBin. See the Tcl configuration
# reference in README.md for parameter documentation.

# Leaflet assignment
set leaflet_sorting_algorithm 2
set leaflet_sorter_2_reference_sel "none"
set leaflet_reassign_interval 5
set restrict_leaflet_sorter_to_Rmax 0

# Coordinate preparation
set center_and_align 0
set use_qwrap 0

# Helix assignment script
set helix_assignment_script "assign_helices_GPCR_general.tcl"

# Frames and trajectory sampling
set start_frame 0
set dt 1

# Protein selections and helix assignment
set backbone_selstr "name BB"
set protein_selstr "name BB SC1 to SC4"
set chainlist [list A]
set helixlist [list 1 2 3 4 5 6 7]
set midplane_selstr "occupancy 1 to 7"

# Lipid selections and output names
set atomsels [list "resname CHOL"]
set filename_stems [list "CHOL"]
set headnames [list "name ROH"]
set tailnames [list "name C2"]

# Polar grid
set Rmax 20.
set Rmin 0.
set dr 1.
set Ntheta 50
