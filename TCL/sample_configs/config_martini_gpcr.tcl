# Configuration example for a GPCR-like Martini system.
#
# Edit the atom selections, paths, chain names, and helix assignments to match
# the loaded molecule before running polarDensityBin.

# Leaflet assignment:
#   0 = compare the specified head and tail selections. headnames and
#       tailnames are atom-selection strings and may select multiple beads.
#   1 = compare the default lipid termini (legacy method)
#   2 = compare lipid position with leaflet_sorter_2_reference_sel
#   3 = compare each lipid with a local membrane midplane
#   4 = assign every lipid to the upper leaflet
set leaflet_sorting_algorithm 2
set leaflet_sorter_2_reference_sel "none"
set leaflet_reassign_interval 5
set restrict_leaflet_sorter_to_Rmax 0

# Coordinate preparation.
# Set center_and_align to 1 only when the protein should be recentered and
# aligned before binning. qwrap is an optional alternative wrapping method;
# it requires an orthorhombic unit cell and is disabled by default.
# A pre-compiled qwrap library is currently included with DTA.
set center_and_align 0
set use_qwrap 0

# Set this to the path of the helix-assignment script you prepared for the
# loaded system. The utilities path is configured in polarDensity_for_DTA.tcl.
set helix_assignment_script "assign_helices_GPCR_general.tcl"

# Frames and trajectory sampling. end_frame defaults to the final loaded frame.
set start_frame 0
set end_frame [expr {[molinfo top get numframes] - 1}]
set dt 1

# Protein selections and helix assignment. The helix-assignment script assigns
# occupancy values to the transmembrane helices; helixlist contains those
# occupancy values, not atom indices. chainlist contains chain identifiers
# used when writing protein coordinates. If all subunits share one chain, use
# a one-item chainlist.
set backbone_selstr "name BB"
set protein_selstr "name BB SC1 to SC4"
set chainlist [list A]
set helixlist [list 1 2 3 4 5 6 7]
set midplane_selstr "occupancy 1 to 7"

# Lipid selections and output names. There must be one filename stem for
# every atom selection. headnames and tailnames are used when algorithm 0 is
# selected. Each is an atom-selection string, and can select multiple beads.
set atomsels [list "resname CHOL"]
set filename_stems [list "CHOL"]
set headnames [list "name ROH"]
set tailnames [list "name C2"]

# Polar grid. (Rmax - Rmin) must be evenly divisible by dr.
set Rmax 20.
set Rmin 0.
set dr 1.
set Ntheta 50
