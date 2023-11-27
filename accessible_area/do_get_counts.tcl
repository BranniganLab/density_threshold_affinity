source get_counts.tcl
set UTILS "../../densitymap/TCL/helpers"
set USE_QWRAP 0
source ../polarDensity_for_DTA.tcl
set ASSIGN_LEAFLETS 1
set species "DPPC"
set area 154

set step [expr $area**(0.5)]

set all [atomselect top "name PO4"]
set XMIN [lindex [measure minmax $all] 0 0]
set YMIN [lindex [measure minmax $all] 0 1]
set XMAX [lindex [measure minmax $all] 1 0]
set YMAX [lindex [measure minmax $all] 1 1]

set XMIN [expr $XMIN+$step]
set YMIN [expr $YMIN+$step]
set XMAX [expr $XMAX-$step]
set YMAX [expr $YMAX-$step]

if {$ASSIGN_LEAFLETS == 1} {
	set nframes [molinfo top get numframes]
	set all_lipids [atomselect top "resname $species"]
	set resids [lsort -unique [$all_lipids get resid]]
	foreach resid $resids {
		for {set frm 0} {$frm < $nframes} {incr frm} {
			local_mid_plane "resname $species and resid $resid" $frm
		}
	}
}

set box_width [lindex [measure minmax $all] 1 0]
$all delete

set outfile [open "counts_${area}_GBsorted.out" w]   

foreach bfield {1 '-1'} {
	set xmin $XMIN
	while {$xmin < $XMAX} {
		set ymin $YMIN
		while {$ymin < $YMAX} {
			puts "Running at $xmin $ymin leaflet $bfield"
			set data [get_count_with_area $area $xmin $ymin "resname DPPC and beta $bfield"]
			puts $outfile $data

			set ymin [expr $ymin + $step]
		}
		set xmin [expr $xmin + $step]
	}
}

close $outfile
