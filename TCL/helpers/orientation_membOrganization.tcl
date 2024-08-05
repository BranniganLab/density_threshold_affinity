
# namespace explicitly dealing with how the membrane is organized/centerd
namespace eval membOrganization {

	;# Centers around the protien and centers the system at 0,0,0
	# sharp changed atsel_str to atsel_str 8/5
	proc Center { atsel_str } {
	    pbc wrap -center com -centersel ${atsel_str} -orthorhombic  -all
	    set nframes [molinfo top get numframes]
	    set zero [veczero]
	    for {set frames 0} {$frames < $nframes} {incr frames} {
	        set sel [atomselect top "all" frame $frames]
	        set com [measure center $sel]
	        $sel delete
	        set sell [atomselect top "all" frame $frames]
	        set mov [vecsub $zero $com]
	        $sell moveby $mov
	        $sell delete
	    }
	}

	;# Alignment based off vmd alignment
	proc Align { atsel_str } {
	    set nframes [molinfo top get numframes]
	    set ref [atomselect top $atsel_str frame 0]
	    for {set frames 1} {$frames < $nframes} {incr frames} {
	        set com [atomselect top $atsel_str frame $frames]
	        set TM [measure fit $com $ref]
	        $com delete
	        set move_sel [atomselect top "all" frame $frames]
	        $move_sel move $TM
	        $move_sel delete
	    }
	    $ref delete
	}

	;#determines the average acyl chain length for a given species 
	;#usage unclear
	proc avg_acyl_chain_len {species acylchain_selstr} {
	    set acyl_num 0
	    set sel [atomselect top "$species"]
	    set sel_resname [lsort -unique [$sel get resname]]
	    #set sel_num [llength [lsort -unique [$sel get resname]]]
	    $sel delete
	    foreach res $sel_resname {
	        set sel [atomselect top "${species} and (resname $res) and $acylchain_selstr"]
	        set sel_len [llength [lsort -unique [$sel get name]]]
	        # 6 is the longest chain in Martini
	        # If there is a chain longer -> the lipid is
	        # a homoacid, need to add another value to 
	        # divide by
	        if {$sel_len > 6} {
	            lappend sel_resname "${res}"
	        }
	        $sel delete
	        set acyl_num [expr $acyl_num + $sel_len]
	    }
	    set avg_acyl_chain [ expr (1.0 * $acyl_num / [llength $sel_resname]) ]
	    if {$avg_acyl_chain < 1} {
	        return 1
	    }
	    return $avg_acyl_chain
	    
	}

	;# confirms your box is either square or paraelleogram-ish
	;# will perform qwrap or pbc wrap depending
	# sharp changed input to atsel_str 8/5
	proc center_and_wrap_system {atsel_str} {
	    global params
	    set pbc_angles [molinfo top get {alpha beta gamma}]
	    
	    set sel [atomselect top "$atsel_str" frame 0]
	    set com [measure center $sel weight mass]
	    
	    set counter_i 0
	    # continues to try and recenter's box until ~ 0'ed out
	    while {[expr abs([lindex $com 0])] > 1.0 &&  [expr abs([lindex $com 1])] > 1.0} {
	        
	        if {$counter_i > 5} {
	            puts "Script was unable to converge system to (0,0,0)"
	            puts "Please check your system vissually, there may be"
	            puts "unintended artifacts"
	            $sel delete
	            return
	        }
	        
	        if {($params(use_qwrap)==0) || (([lindex $pbc_angles 0]!=90.0) && ([lindex $pbc_angles 0]!=90.0) && ([lindex $pbc_angles 0]!=90.0))} {
	            puts "qwrap may not be optimal for your system...\n"
	            puts "Running pbc wrap. To verify proper centering"
	            puts "pbc wrap will be run multiple times" ; after 100
	            foreach i {0 1 2 3} {
	                pbc wrap -centersel "$atsel_str" -all
	            }
	        } else {
	            qwrap centersel "$atsel_str" ;#center entire system at ~0,0,0
	        }
	        set com [measure center $sel weight mass]
	        incr counter_i
	    }
	    $sel delete
	}
}

# namespace dealing explicitly with lipid assessment
# which sorter, assigning lipid to leaflets...
namespace eval lipid_assessment {
	;#determines leaflet based on relative height of specified head and tail beads
	proc leaflet_sorter_0 {atsel_in head tail frame_i} {
	    
	    #puts "Sorting into leaflets using leaflet_sorter_0"
	    set sel_resid [atomselect top "$atsel_in" frame $frame_i]
	    set sel_head [atomselect top "$atsel_in and name $head" frame $frame_i]
	    set sel_tail [atomselect top "$atsel_in and name $tail" frame $frame_i]
	    
	    set head_Z [${sel_head} get z] 
	    set tail_Z [${sel_tail} get z] 
	    # set head_Z [expr abs( [lindex ${head_Z} 0] )]
	    # set tail_Z [expr abs( [lindex ${tail_Z} 0] )]

	    if {$head_Z < $tail_Z } { 
	        $sel_resid set user2 -1
	        return -1 
	    } else { 
	        $sel_resid set user2 1
	        return 1 
	    }
	    $sel_resid delete
	    $sel_head delete
	    $sel_tail delete
	}

	;#originally by Liam Sharp; procedure that was used in JCP 2021 for nAChR
	;#similar to leaflet_sorter_0 but autoselects head and tail beads 
	proc leaflet_sorter_1 {atsel_in frame_i} {
	    #puts "Sorting into leaflets using leaflet_sorter_1"
	    set sel_resid [atomselect top "$atsel_in" frame $frame_i]
	    set ind 1
	    if { [string range [lsort -unique [$sel_resid get resname]] end-1 end] == "PA" } {
	        set ind 0
	    }
	    set sel_Z [${sel_resid} get z] 
	    if {[lindex ${sel_Z} $ind] < [lindex ${sel_Z} end] } { 
	        $sel_resid set user2 -1
	        return -1 
	    } else { 
	        $sel_resid set user2 1
	        return 1 
	    }
	    $sel_resid delete
	}

	;#originally by Jahmal Ennis, designed for cholesterol 
	proc leaflet_sorter_2 {atsel_in frame_i} {

	    #puts "Sorting into leaflets using leaflet_sorter_2"
	    set sel_resid [atomselect top "$atsel_in" frame $frame_i]
	    set ind 1
	    if { [string range [lsort -unique [$sel_resid get resname]] end-1 end] == "PA" } {
	        set ind 0
	    }
	    if { [lsort -unique [$sel_resid get resname]] == "CHOL" } {
	        set ind 2
	    }
	    if { $ind == 2 } {
	        set chol_com_z [lindex [measure center $sel_resid weight mass] end]
	        if {$chol_com_z < 0} {
	            $sel_resid set user2 -1
	            return -1
	        } else {
	            $sel_resid set user2 1
	            return 1
	        }
	    }
	}


	;# Determines if the lipid is in the outer or inner leaflet and sets the user2 value accordingly
	;# Algorithm is determined by user: 
	;# 0: determines leaflet based on relative height of specified head and tail beads
	;# 1: originally by Liam Sharp; procedure that was used in JCP 2021 for nAChR; similar to leaflet_sorter_0 but autoselects head and tail beads; more appropriate for situations with many species
	;# 2: originally by Jahmal Ennis, determines whether the auto-determined headbead is above or below the center of mass (of what? the system?); more appropriate for rigid lipids like cholesterol that frequently invert or lie at parallel to the membrane
	proc leaflet_detector {atsel_in head tail frame_i leaflet_sorting_algorithm} {
	    if {$leaflet_sorting_algorithm == 0} {
	        leaflet_sorter_0 $atsel_in $head $tail $frame_i
	    } elseif { $leaflet_sorting_algorithm == 1 } {
	        leaflet_sorter_1 $atsel_in $frame_i
	    } elseif { $leaflet_sorting_algorithm == 2 } {
	        leaflet_sorter_2 $atsel_in $frame_i
	    } else { 
	        puts "Option $leaflet_sorting_algorithm not recognized as a leaflet sorting option.  Defaulting to option 1." 
	    }
	}


	;# Calculates the total number of lipids and beads of the given species in each leaflet 
	;# Assigns the leaflet to user2 
	;# Returns the following list : [["lower" lower_leaflet_beads lower_leaflet_lipids] ["upper" upper_leaflet_beads upper_leaflet_lipids]] 
	proc frame_leaflet_assignment {species headname tailname lipidbeads_selstr frame_i frame_f} {
	    global params
	    set sel [ atomselect top "(($species)) and $lipidbeads_selstr"  frame $frame_i]
	    set sel_num [llength [lsort -unique [$sel get resid] ] ]
	    set sel_resid_list [lsort -unique [$sel get resid] ]
	    set totals {}
	    if {$sel_num < 1} {
	        set totals [[list 0 0] [list 0 0]] 
	    } else {
	        #assign leaflets from $frame_i to user2 field of each bead for this species
	        foreach sel_resid $sel_resid_list {
	            set selstring "${species} and (resid $sel_resid) and $lipidbeads_selstr"
	            set leaflet [leaflet_detector $selstring $headname $tailname $frame_i $params(leaflet_sorting_algorithm)]
	        }
	        #copy leaflet values from $frame_i to all frames between $frame_i and $frame_f
	        set leaflet_list [$sel get user2] 
	        for {set interim_frame [expr $frame_i + 1]} {$interim_frame < [expr $frame_f]} {incr interim_frame} {
	            $sel frame $interim_frame
	            $sel set user2 $leaflet_list
	        }
	        #count the number of lipids and the number of beads in each leaflet
	        foreach leaf [list  "(user2<0)" "(user2>0)"] txtstr [list "lower" "upper"] {
	            set leaf_sel [ atomselect top "(${species} and $leaf)"  frame $frame_i]
	            set num_beads [$leaf_sel num]
	            set num_lipids [llength [lsort -unique [$leaf_sel get resid] ]]
	            lappend totals [list $txtstr $num_beads $num_lipids]
	            $leaf_sel delete
	        }
	    }
	    $sel delete
	    return $totals
	}

	;# Calculates the total number of lipids and beads of the given species in each leaflet 
	;# Returns the following list : [["lower" lower_leaflet_beads lower_leaflet_lipids] ["upper" upper_leaflet_beads upper_leaflet_lipids]] 
	proc trajectory_leaflet_assignment {species headname tailname lipidbeads_selstr} { 
	    global params
	    set num_reassignments 0
	    for {set update_frame $params(start_frame)} {$update_frame < $params(end_frame)} {incr update_frame $params(dt)} {
	        frame_leaflet_assignment $species $headname $tailname $lipidbeads_selstr $update_frame [expr $update_frame + $params(dt)] 
	        incr num_reassignments
	    }
	    puts "Checked for leaflet reassignments $num_reassignments times."
	}
	    
	;#Reinitializes the user2 value for selected beads in selected frames 
	proc clean_leaflet_assignments {species lipidbeads_selstr} {
	    global params
	    set sel [ atomselect top "$species and $lipidbeads_selstr"]
	    set selnum [$sel num]

	    for {set update_frame $params(start_frame)} {$update_frame < ${ _frame}} {incr update_frame} {
	        $sel frame $update_frame
	        $sel set user2 [lrepeat $selnum 0.0]
	        puts "Cleaning $selnum beads of leaflet assignments in frame $update_frame"
	    }
	    $sel delete
	}

}