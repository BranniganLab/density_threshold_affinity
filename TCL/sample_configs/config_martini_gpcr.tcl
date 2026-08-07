set currentPath [file normalize [info script]]

# Options for leaflet_sorting_algorithm
# 0 for sorting based on orientation of specified head & tail; 
# 1 for legacy sorting based on orientation of default termini (aka classic local_mid_plane); 
# 2 for sorting based on position relative to origin (tested only with cholesterol so far) 
set leaflet_sorting_algorithm 2; 
set restrict_leaflet_sorter_to_Rmax 1;

set center_and_align 1
set use_qwrap 0
set utils "${currentPath}/../utilities" 

set dt 1
set leaflet_reassign_interval 1; #how frequently to reassign lipids to leaflets
#set start_frame 0 ; #optional
#set end_frame 10  ; #optional 

set backbone_selstr "name BB" ;#selection string used to define the protein backbone
set protein_selstr "name BB SC1 to SC4" ;#selection string used to define the entire protein

set atomsels [list "resname CHOL and name ROH"]
set filename_stems [list "CHOL_ROH"]

set chainlist [list A] ;#list of chain names for the protein
set helixlist [list 1 2 3 4 5 6 7]; #indices for individual secondary structure elements 
set helix_assignment_script "${currentPath}/../sample_helix_assignment_scripts/assign_helices_GPCR_general.tcl" ;# script that will assigns occupancies in helixlist to different secondary structure elements 
set midplane_selstr "occupancy 1 to 7" ;# selection that includes all transmembrane helices

set Rmax 20. ;##maximum radius of polar density map
set Rmin 0. ;#minimum radius
set dr 1 ;#radial bin width 
set Ntheta 50; #number of angular bins 

