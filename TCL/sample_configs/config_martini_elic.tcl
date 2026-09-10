# Configuration example based on a Martini ELIC system.
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
set leaflet_sorting_algorithm 1
set leaflet_sorter_2_reference_sel "none"
# leaflet_reassign_interval controls how often leaflet identities are
# recalculated. Between reassignment frames, the previous identities are reused.
set leaflet_reassign_interval 1
set restrict_leaflet_sorter_to_Rmax 0

# Coordinate preparation.
# Set center_and_align to 1 only when the protein should be recentered and
# aligned before binning. qwrap is an optional alternative wrapping method;
# it requires an orthorhombic unit cell and is disabled by default.
# A pre-compiled qwrap library is bundled with DTA and is loaded automatically
# when use_qwrap is enabled.
set center_and_align 0
set use_qwrap 0

# Set this to the path of the helix-assignment script you prepared for the
# loaded system.
set helix_assignment_script "assign_helices_ELIC_general.tcl"

# Frames and trajectory sampling. end_frame defaults to the final loaded frame.
set start_frame 0
# Optional: set end_frame to a specific final frame. If omitted, the final
# loaded frame is used.
# dt controls the interval, in frames, between frames included in the density
# calculation. dt=1 analyzes every frame.
set dt 1

# Protein selections and helix assignment. The helix-assignment script assigns
# occupancy values to the transmembrane helices; helixlist contains those
# occupancy values, not atom indices. chainlist contains the chain identifiers
# used to distinguish protein subunits/components when writing inclusion
# coordinates. Each subunit that should be represented separately must have a
# distinct chain ID.
set backbone_selstr "name BB"
set protein_selstr "name BB SC1 to SC4"
set chainlist [list A B C D E]
set helixlist [list 1 2 3 4]
# Selection used to define the membrane midplane when protein/inclusion
# coordinates are divided between the upper and lower leaflet outputs. Its
# mass-weighted z center is averaged over the analyzed trajectory.
set midplane_selstr "occupancy 1 to 4"

# Lipid selections and output names. There must be one filename stem for
# every atom selection. headnames and tailnames are used only by
# leaflet_sorting_algorithm 0. Each entry is an atom-selection string for one
# lipid species. If a selection contains multiple beads, their z coordinates
# are averaged before the mean head and tail heights are compared.
set atomsels [list "resname POPG"]
set filename_stems [list "POPG"]
set headnames [list "name PO4"]
set tailnames [list "name C4"]

# Polar grid. (Rmax - Rmin) must be evenly divisible by dr.
set Rmax 20.
set Rmin 0.
set dr 5.
set Ntheta 50
