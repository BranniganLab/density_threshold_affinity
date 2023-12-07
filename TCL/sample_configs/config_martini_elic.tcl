# Options for LEAFLET_SORTING_ALGORITHM
# 0 for sorting based on orientation of specified head & tail; 
# 1 for legacy sorting based on orientation of default termini (aka classic local_mid_plane); 
# 2 for sorting based on position relative to origin (tested only with cholesterol so far) 
set LEAFLET_SORTING_ALGORITHM 1; 
set CENTER_AND_ALIGN 0
set USE_QWRAP 0
set UTILS "./helpers" 

set dt 1
set leaflet_reassign_t 5; #how frequently to reassign lipids to leaflets
#set start_frame 0 ; #optional
#set end_frame 10  ; #optional 

set backbone_selstr "name BB" ;#selection string used to define the protein backbone
set protein_selstr "name BB SC1 to SC4" ;#selection string used to define the entire protein

set lipids [list "POPG"] ;# list of all species to bin
set headnames [list "PO4"] ; #lists one headgroup atom/bead name per lipid species 
set tailnames [list "C4"] ; #lists one terminal atom/bead name per lipid species 
set lipidbeads_selstrs [list "all"]; #lists one selection string per lipid species; indicates which lipid atom/beads should be counted in the density plot 
set acylchain_selstrs $lipidbeads_selstrs; #list of beads used to determine chain length; same format as lipidbeads_selstrs. Was originally intented to only hold selection strings containing the acyl chain beads.  

set chainlist [list A B C D E] ;#list of chain names for the protein
set helixlist [list 1 2 3 4]; #indices for individual secondary structure elements 
set helix_assignment_script assign_helices_ELIC_general.tcl ;# script that will assigns occupancies in helixlist to different secondary structure elements 
set midplane_selstr "occupancy 1 to 4" ;# selection that includes all transmembrane helices

set Rmax 20. ;##maximum radius of polar density map
set Rmin 0. ;#minimum radius
set dr 5 ;#radial bin width 
set Ntheta 50; #number of angular bins 

