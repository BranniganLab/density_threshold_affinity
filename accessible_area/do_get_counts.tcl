source get_counts.tcl
set xmin 10
set ymin 10

set area 154
set data [get_count_with_area $area $xmin $ymin "resname DPPC"]

set outfile [open "counts_${area}.out" w+]    
puts $outfile $data
close $outfile