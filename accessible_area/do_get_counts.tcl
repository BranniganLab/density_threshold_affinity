variable help_me_dir [file dirname [file normalize [info script]]]
source $help_me_dir/utilities/TCL/get_counts.tcl
source $help_me_dir/../TCL/polarDensity_for_DTA.tcl


proc run_leaflet_sorter {atsel sorter_num {stride 1}} {
	set all_lipids [atomselect top $atsel]
	set resids [lsort -integer -unique [$all_lipids get resid]]
	set nframes [molinfo top get numframes]
	$all_lipids delete
	foreach resid $resids {
		if {[expr $resid%100] == 0} { puts "on resid $resid at [clock seconds]" }
		set seltext "$atsel and resid $resid"
		for {set frm 0} {$frm < $nframes} {incr frm $stride} {
			if {$sorter_num == 1} {
				leaflet_sorter_1 $seltext $frm	
			} elseif {$sorter_num == 3} {
				leaflet_sorter_3 $seltext $frm	
			}
			
		}
	}
}


set atsels [list "resname CHOL"]
set names [list "CHOL"]

set ASSIGN_LEAFLETS 1
set leaflet_sorter_option 1
set GET_COUNTS 1

set area 85
set stride 10
set fieldid "user2"


if {$ASSIGN_LEAFLETS == 1} {
	if {[lsearch -exact "1 3" $leaflet_sorter_option] == -1} {
        	puts "Option ${leaflet_sorter_option} not recognized as a leaflet sorting option. Defaulting to option 1."
		set leaflet_sorter_option 1
	}
	foreach atsel $atsels {
		puts "assigning leaflets for $atsel with leaflet sorter $leaflet_sorter_option"
		run_leaflet_sorter $atsel $leaflet_sorter_option $stride
	}
}


if {$GET_COUNTS == 1} {
	set all [atomselect top "name PO4"]
	set minmax [measure minmax $all]
	$all delete

	set XMIN [lindex $minmax 0 0]
	set YMIN [lindex $minmax 0 1]
	set XMAX [lindex $minmax 1 0]
	set YMAX [lindex $minmax 1 1]

	set step [expr $area**(0.5)]
	set XMIN [expr $XMIN+$step]
	set YMIN [expr $YMIN+$step]
	set XMAX [expr $XMAX-$step]
	set YMAX [expr $YMAX-$step]

	foreach atsel $atsels name $names {
		set outfile [open "${name}_counts_${area}.out" w]   
		foreach leaflet_id {1 '-1'} {
			set xmin $XMIN
			while {$xmin < $XMAX} {
				set ymin $YMIN
				while {$ymin < $YMAX} {
					puts "Running at $xmin $ymin leaflet $leaflet_id"
					set data [get_count_with_area $area $xmin $ymin "(${atsel}) and $fieldid $leaflet_id" top 0 -1 $stride]
					puts $outfile $data
					set ymin [expr $ymin + $step]
				}
				set xmin [expr $xmin + $step]
			}
		}
		close $outfile
	}
}
