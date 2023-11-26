
set CENTER_AND_ALIGN 0
set USE_QWRAP 0

set UTILS "../densitymap/TCL/helpers" 
source ./polarDensity_for_DTA.tcl

set HeadNames [atomselect top "name PO4"] ;#B2 is DDM
set lipids [lsort -unique [$HeadNames get resname]]
puts $lipids
$HeadNames delete
puts $lipids	
set Rmax 60.
set Rmin 0.
set dr 5.
set Ntheta 10
set dt 1
set sample_frame 0

source ${UTILS}/assign_helices_ELIC_CG.tcl
foreach lip $lipids {
	polarDensityBin $lip.dat $lip $Rmin $Rmax $dr $Ntheta $dt $sample_frame	[list A B C D E] [list 1 2 3 4]
}
