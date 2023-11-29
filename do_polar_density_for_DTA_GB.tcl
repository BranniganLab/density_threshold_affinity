
set CENTER_AND_ALIGN 0
set USE_QWRAP 0

set UTILS "../densitymap/TCL/helpers" 
source ./polarDensity_for_DTA.tcl
set Rmax 35.
set Rmin 0.
set dr 5
set Ntheta 10
set dt 1
set sample_frame 0



#martini
#set midplane_selstr "name PO4 ROH C3 PO41" ;# used in z_mid
#set backbone_selstr "name BB" ;#used in Protein_Position
#set acylchain_selstr "(not name NH3 NC3 GL1 GL2 AM1 AM2 PO4 CNO CN0 C1 C2 C3)" ;#used in acyl_chain_length
#set lipidbeads_selstr "(not name PO4)"; #selection for beads to use in the density calculation for species $species
#set protein_selstr "name BB SC1 to SC4" ;#martini
#set HeadNames [atomselect top "name ROH"] ;#B2 is DDM

#charmm
set midplane_selstr "name P" ;# used in z_mid
set backbone_selstr "alpha" ;#used in Protein_Position
set protein_selstr "protein" ; #used in helix assignment
set acylchain_selstr "(not name NH3 NC3 GL1 GL2 AM1 AM2 PO4 CNO CN0)" ;#used in acyl_chain_length, should select every atom
set lipidbeads_selstr "(not name PO4)"; #selection for beads to use in the density calculation for species $species; should select every atom
set HeadNames [atomselect top "name P"] ;#B2 is DDM

#GPCR settings
#set helix_assignment_script "assign_helices_GPCR_general.tcl";#gpcr 
#set chainlist [list A];#gpcr
#set helixlist [list 1 2 3 4 5 6 7];#gpcr

#pLGIC settings
set chainlist [list A B C D E] ;#plgic
set helixlist [list 1 2 3 4];#plgic
set helix_assignment_script assign_helices_ELIC_general.tcl

set lipids [list "POPG"]	

source ${helix_assignment_script}
foreach lip $lipids {
	polarDensityBin $lip.dat $lip $Rmin $Rmax $dr $Ntheta $dt $sample_frame	$chainlist $helixlist $midplane_selstr $backbone_selstr $acylchain_selstr $lipidbeads_selstr
}

