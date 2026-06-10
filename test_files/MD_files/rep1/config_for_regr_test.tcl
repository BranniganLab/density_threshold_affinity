# Options for leaflet_sorting_algorithm:
# 0 for sorting based on orientation of specified head & tail; 
# 1 for legacy sorting based on orientation of default termini (aka classic local_mid_plane); 
# 2 for sorting based on position relative to origin (tested only with cholesterol so far) 
set leaflet_sorting_algorithm 3; 
set restrict_leaflet_sorter_to_Rmax 1;

set center_and_align 1
set use_qwrap 0
set utils "../../../../TCL/utilities" 

set dt 1
#set leaflet_reassign_interval 1; #optional
#set start_frame 0 ; #optional
#set end_frame 10  ; #optional 

set backbone_selstr "name BB" ;#selection string used to define the protein backbone
set protein_selstr "name BB SC1 to SC4" ;#selection string used to define the entire protein

set atomsels [list "resname CDL1"] ;# list of all species to bin
set filename_stems [list "regr_test"]

set chainlist [list A B C D E] ;#list of chain names for the protein
set helixlist [list 1 2 3 4]; #indices for individual secondary structure elements 
set helix_assignment_script "../assign_helices_ELIC_CG.tcl" ;# script that will assigns occupancies in helixlist to different secondary structure elements 
set midplane_selstr "occupancy 1 to 4" ;# selection that includes all transmembrane helices

set Rmax 50. ;##maximum radius of polar density map
set Rmin 0. ;#minimum radius
set dr 5 ;#radial bin width 
set Ntheta 90; #number of angular bins 
