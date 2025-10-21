proc get_counts { xmin ymin xmax ymax seltext {molid top} {start 0} {end -1} {stride 1}} {
	set sel [atomselect $molid "$seltext and x > $xmin and x < $xmax and y > $ymin and y < $ymax"]
	if {$end > 0} {
		set stop $end
	} else {
		set stop [molinfo top get numframes]
	}
	set data ""
	for {set frm $start} {$frm < $stop} {incr frm $stride} {
		$sel frame $frm 
		$sel update
		set toadd [$sel num]
		lappend data $toadd
	}
	$sel delete
	return $data
}

proc get_count_with_area { area ymin xmin seltext {molid top} {start 0} {end -1} {stride 1} } {
	set width [expr "$area**(0.5)"]
	set xmax [expr $xmin+$width]
	set ymax [expr $ymin+$width]

	set data [get_counts $xmin $ymin $xmax $ymax $seltext $molid $start $end $stride]
	return $data
}
