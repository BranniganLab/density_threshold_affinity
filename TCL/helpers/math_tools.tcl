# name space for math function used not default in TCL

namespace eval mathTools {
# get_avg_area
	#
	# Calculates the average area of a molecule in a simulation box.
	# Arguments:
	#     molid (str): The molecule ID of the molecule for which the average area needs to be calculated.
	# Results:
	#     float: The average area of the molecule in the simulation box.
	proc get_avg_area {molid} {
	    set box [pbc get -all]
	    set xbox [list]
	    set ybox [list]
	    foreach r $box {
	        lappend xbox [lindex $r 0]
	        lappend ybox [lindex $r 1]
	    }
	    set x [expr 1.0 * [vecsum $xbox] / [llength $xbox]]
	    set y [expr 1.0 * [vecsum $ybox] / [llength $ybox]]
	    set avg [expr 1.0 * $x * $y]
	    return $avg
	}
	# RtoD
	#
	# Converts an angle from radians to degrees
	# Arguments:
	#   float: an angle in radians
	# Outputs:
	#   float: the same angle in degrees
	proc RtoD {r} {
	    set pi 3.14159265358979323846
	    return [expr $r*180.0/$pi]
	}

	#Converts degrees to radians
	proc DtoR {d} {
	    set pi 3.1415926535897931
	    set numerator [expr $d*$pi]
	    set out [expr $numerator/180.0]
	    return $out
	}
				

	# get_theta
	#
	# Gets the angle from the (x,y) coordinates
	# Arguments:
	#   float: x position
	#   float: y position
	# Outputs:
	#   float: the angle of the point in degrees relative to the origin
	proc get_theta {x y} {
	    set pi 3.14159265358979323846
	    set tmp  [expr {atan2($y,$x)}]
	    if {$tmp < 0} {
	        set theta [expr 2*$pi + $tmp]    
	    } else {
	        set theta $tmp
	    }
	    return [RtoD $theta]
	}
}