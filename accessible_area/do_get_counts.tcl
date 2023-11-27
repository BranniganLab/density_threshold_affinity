source get_counts.tcl

set species "DPPC"
set midresid 1522
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

set upper [atomselect top "resname $species and resid < $midresid"]
$upper set beta 1
set lower [atomselect top "resname $species and resid >= $midresid"]
$lower set beta -1


set box_width [lindex [measure minmax $all] 1 0]
$all delete

set outfile [open "counts_${area}.out" w]   

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
