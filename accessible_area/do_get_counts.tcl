variable help_me_dir [file dirname [file normalize [info script]]]
source $help_me_dir/get_counts.tcl
set UTILS "$help_me_dir/../TCL/helpers"
set USE_QWRAP 0
source $help_me_dir/../TCL/polarDensity_for_DTA.tcl

proc use_old_sorter {species {stride 1}} {
	set all_lipids [atomselect top "resname $species"]
	set resids [lsort -integer -unique [$all_lipids get resid]]
	set nframes [molinfo top get numframes]
	foreach resid $resids {
		if {[expr $resid%100] == 0} { puts "on resid $resid at [clock seconds]" }
		set seltext "resname $species and resid $resid"
		for {set frm 0} {$frm < $nframes} {incr frm $stride} {
			leaflet_sorter_3 $seltext $frm
		}
	}
}

proc assign_all_frames {species {stride 1}} {
	set nframes [molinfo top get numframes]
	set all_lipids [atomselect top "resname $species"]
	set resids [lsort -integer -unique [$all_lipids get resid]]
	$all_lipids delete
	foreach resid $resids {
		if {[expr $resid%100] == 0} { puts "on resid $resid at [clock seconds]" }
		for {set frm 0} {$frm < $nframes} {incr frm $stride} {
			leaflet_sorter_1 "resname $species and resid $resid" $frm
		}
	}
}


set ASSIGN_LEAFLETS 1
set GET_COUNTS 1
set spp [list "DPPC"]
set area 52
set stride 10
set fieldid "user2"

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
	set fieldid "user2"
	foreach species $spp {
		puts "assigning leaflets of $species"
		#use_old_sorter $species $stride
		assign_all_frames $species $stride
	}
}

set box_width [lindex [measure minmax $all] 1 0]
$all delete

foreach species $spp {
	if {$GET_COUNTS == 1} {
		set outfile [open "${species}_counts_${area}.out" w]   

		foreach field {1 '-1'} {
			set xmin $XMIN
			while {$xmin < $XMAX} {
				set ymin $YMIN
				while {$ymin < $YMAX} {
					puts "Running at $xmin $ymin leaflet $field"
					set data [get_count_with_area $area $xmin $ymin "resname $species and $fieldid $field" top 0 -1 $stride]
					puts $outfile $data

					set ymin [expr $ymin + $step]
				}
				set xmin [expr $xmin + $step]
			}
		}

		close $outfile
	}
}
