source get_counts.tcl
set xmin 0

set area 154

set step [expr $area**(0.5)]
set box_width [lindex [pbc get] 0 0]
set usable_width [expr $box_width - $step]

set outfile [open "counts_${area}.out" w+]   


while {$xmin < $usable_width} {
	set ymin 0
	while {$ymin < $usable_width} {
		puts "Running at $xmin $ymin"
		set data [get_count_with_area $area $xmin $ymin "resname DPPC"]
		puts $outfile $data

		set ymin [expr $ymin + $step]
	}

	set xmin [expr $xmin + $step]
}

close $outfile