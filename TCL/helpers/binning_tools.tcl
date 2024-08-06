# namespace for tools that assist with binning... 
# a few might need to be moved...

namespace eval binningTools {

	
	# lcount 
	#
	# Takes a list as input and returns a list of unique elements in the input list along with their counts.
	# Arguments:
	#   list: A list of elements.
	# Outputs:
	#   A list containing two lists:
	#   - The first list contains the unique elements from the input list.
	#   - The second list contains the counts of each unique element in the input list.
	proc lcount list {
	    foreach x $list {lappend arr($x) {}}
	    set res1 {}	
	    set res2 {}
	    foreach name [array names arr] {
	        lappend res1 $name
	        lappend res2 [llength $arr($name)]
	    }
	    set res [list $res1 $res2]
	    return $res
	}


	proc  Coord_Profile {Xi Xf Yi Yf} {
	    if {($Xi < $Xf) && ($Yi < $Yf)} {
	        set a "x > $Xi and x < $Xf and y > $Yi and y < $Yf"
	    } elseif { ($Xi > $Xf) && ($Yi < $Yf)} {
	        set a "x < $Xi and x > $Xf and y > $Yi and y < $Yf"
	    } elseif {($Xi < $Xf) && ($Yi > $Yf)} {
	        set a "x > $Xi and x < $Xf and y < $Yi and y > $Yf"
	    } elseif {($Xi > $Xf) && ($Yi > $Yf)} {
	        set a "x < $Xi and x > $Xf and y < $Yi and y > $Yf"
	    }
	    return $a
	}

	proc Time_Predict {step_list estimated_steps} { ;#linear, is not really accurate
	    set avg 0
	    foreach v $step_list { set avg [expr $avg + $v] }
	    set avg [expr 1.0*$avg/[llength $step_list]]
	    set predict [expr 1.0*$avg*$estimated_steps/60.0]
	    puts "\n\nEstimated Time to Run: $predict (min)\n\n"
	}

	# z_mid
	#
	# Finds the average mid-plane
	# Arguments:
	#   int: first frame
	#   nframes: The number of frames over which to average
	# Outputs:
	#   float: The average z value of all the beads
	#

	proc z_mid {init_frm nframes} {
	    global params
	    set z_list {}
	    for {set frm ${init_frm}} {${frm} < ${nframes}} {incr frm} {
	        set mid [atomselect top $params(midplane_selstr) frame $frm]
	        lappend z_list [lindex [measure center $mid weight mass] 2]
	        $mid delete
	    }
	    return [expr 1.0*[vecsum $z_list]/([llength $z_list]) ]
	}

	#write radial and theta bin output to file 
	proc output_bins {fl  ri rf bins} {
	    global params
	    puts -nonewline $fl "[format {%0.2f} $ri]  [format {%0.2f} $rf] [format {%0.2f} $params(dtheta)]  " 
	    puts $fl "$bins" 
	}

	#
	proc theta_histogram {singleFrame_lower singleFrame_upper } {
	    global params
	    set theta_bin_out [list]
	    foreach ud [list $singleFrame_lower $singleFrame_upper ] {
	        #cleanup and output 
	        set theta_bin_counts [lcount $ud]
	        #Shell_Test $shel_count $theta_bin_counts
	        set theta_bins {}
	        for {set ti 0} { $ti<=$params(Ntheta)} {incr ti 1} {
	            set tindex [lsearch [lindex $theta_bin_counts 0]  $ti]
	            if { $tindex >= 0} {
	                set frame_count [expr 1.0 * [lindex [lindex $theta_bin_counts 1] $tindex]] 
	            } else { 
	                set frame_count 0.0
	            }
	            lappend theta_bins $frame_count
	        }
	        lappend theta_bin_out $theta_bins
	    }
	    return $theta_bin_out
	}


}