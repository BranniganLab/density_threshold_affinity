
package require pbctools


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

# RtoD
#
# Converts an angle from radians to degrees
# Arguments:
#   float: an angle in radians
# Outputs:
#   float: the same angle in degrees
proc RtoD {r} {
    global M_PI
    return [expr $r*180.0/$M_PI]
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
    global M_PI
    set tmp  [expr {atan2($y,$x)}]
    if {$tmp < 0} {
        set theta [expr 2*$M_PI + $tmp]    
    } else {
        set theta $tmp
    }
    return [RtoD $theta]
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



;#Outputs radial position of TMD centers to file for each leaflet
;#Calculated using the midplane_selstr specified in the config to determine leaflet.
proc output_inclusion_centers {{a ""} } {
    global params
    ;# list for the chain names
    ;# finds the center of the membranes
    set midplane_height [z_mid $params(start_frame) $params(end_frame)]
    ;# calculates the center of mass for subunit alpha helices in both leaflets
    puts "Writing coordinates for [llength $params(chainlist)] chains and [llength $params(helixlist)] helices per chain"
    foreach eq {"<" ">"} eqtxt {"lwr" "upr"} {
        set fout [open "./Protein${a}_coords_${eqtxt}.dat" w]
        puts $fout  "# chain/occupancy: r_val theta_val ... "
        set warning_text ""
        foreach chnm $params(chainlist) {
            foreach occ $params(helixlist) {
                set sel [atomselect top "(chain ${chnm}) and (occupancy $occ and $params(backbone_selstr)) and (z ${eq} $midplane_height)" frame 0]
                if {[$sel num] != 0} {
                    set com [measure center $sel weight mass]
                    $sel delete
                    set x [lindex $com 0]
                    set y [lindex $com 1]
                    set r [expr sqrt($x*$x+$y*$y)]
                    set theta [get_theta $x $y]
                    puts -nonewline $fout "${chnm}/${occ}: $r $theta "
                } else {
                    set warning_text "${warning_text} ${chnm}/${occ}"
                }
            }
            puts $fout ""
        }
        close $fout
        if {$warning_text ne ""} {
            puts "NOTE: the following chain/occupancy pairs were not found in the ${eqtxt} leaflet:"
            puts $warning_text
        }
    }
}

;# will perform qwrap or pbc wrap, depending on settings and box angles
proc center_and_wrap_system {inpt} {
    global params
    set counter_i 0

    # determine if cell is orthorhombic (requirement of qwrap)
    set pbc_angles [molinfo top get {alpha beta gamma}]
    if {(([lindex $pbc_angles 0]==90.0) && ([lindex $pbc_angles 0]==90.0) && ([lindex $pbc_angles 0]==90.0))} {
        set orthorhombic 1
    } else {
        set orthorhombic 0
    }
    
    set sel [atomselect top "$inpt" frame 0]
    set com [measure center $sel weight mass]

    # tries to recenter box until COM is 0'ed out or proc times out (5 tries)
    while {[expr abs([lindex $com 0])] > 1.0 &&  [expr abs([lindex $com 1])] > 1.0} {
        
        ;# time-out trigger
        if {$counter_i > 5} {
            puts "Script was unable to converge system to (0,0,0)"
            puts "Please check your system visually, there may be"
            puts "unintended artifacts"
            $sel delete
            return
        }
        
        if {($params(use_qwrap)==0) || ($orthorhombic!=1) } {
            ;# use pbc wrap and center at origin
            if {$params(use_qwrap)!=0} {
                puts "qwrap requires orthorhombic cells.\n"
                puts "Running pbc wrap instead and centering system at origin."
            } else {
                puts "Running pbc wrap and centering system at origin."
            }
            pbc wrap -center com -centersel ${inpt} -all
            center_at_origin
        } else {
            ;# use qwrap (automatically centers at origin)
            qwrap centersel "$inpt"
        }

        ;# update the COM and counter
        $sel update
        set com [measure center $sel weight mass]
        incr counter_i
    }

    $sel delete
}

;#determines leaflet based on relative height of specified head and tail beads
proc leaflet_sorter_0 {atsel_in head tail frame_i} {
    #puts "Sorting into leaflets using leaflet_sorter_0"
    set sel_resid [atomselect top "$atsel_in" frame $frame_i]
    set resnames [$sel_resid get resname]
    set resid [$sel_resid get resid]
    if {[llength [lsort -unique $resid]] != 1} {
        error "More than one resid returned from atomselection [$sel_resid text]"
    }

    ;# necessary in case user has selected only part of a lipid with $atsel_in
    set sel_head [atomselect top "resname $resnames and resid $resid and $head" frame $frame_i]
    set sel_tail [atomselect top "resname $resnames and resid $resid and $tail" frame $frame_i]
    
    set head_Z [vecmean [${sel_head} get z]]
    set tail_Z [vecmean [${sel_tail} get z]]
    $sel_head delete
    $sel_tail delete
    
    if {$head_Z < $tail_Z } { 
        $sel_resid set user2 -1
        $sel_resid delete
        return -1 
    } else { 
        $sel_resid set user2 1
        $sel_resid delete
        return 1 
    }
}

;#originally by Liam Sharp; procedure that was used in JCP 2021 for nAChR
;#similar to leaflet_sorter_0 but autoselects head and tail beads 
proc leaflet_sorter_1 {atsel_in frame_i} {
    #puts "Sorting into leaflets using leaflet_sorter_1"
    set sel_resid [atomselect top "$atsel_in" frame $frame_i]

    set resnames [$sel_resid get resname]
    set resid [$sel_resid get resid]
    if {[llength [lsort -unique $resid]] != 1} {
        error "More than one resid returned from atomselection [$sel_resid text]"
    }

    ;# necessary in case user has selected only part of a lipid with $atsel_in
    set sel [atomselect top "resname $resnames and resid $resid" frame $frame_i]
    
    ;# if the lipid has a PA headgroup, its 0th (rather than 1st) bead is PO4
    if { [string range [lsort -unique [$sel get resname]] end-1 end] == "PA" } {
        set ind 0
    } else {
        set ind 1
    }

    set sel_Z [$sel get z]
    $sel delete

    if {[lindex $sel_Z $ind] < [lindex $sel_Z end] } { 
        $sel_resid set user2 -1
        $sel_resid delete
        return -1 
    } else { 
        $sel_resid set user2 1
        $sel_resid delete
        return 1 
    }
    
}

;#originally by Jahmal Ennis, designed for cholesterol 
proc leaflet_sorter_2 {atsel_in refsel_in frame_i} { 
    if {$refsel_in eq "none"} {
        set refsel_com_z 0
    } else {
        set refsel [atomselect top "$refsel_in" frame $frame_i]
        set refsel_com_z [lindex [measure center $refsel weight mass] 2]
        $refsel delete
    }
    set lipidsel [atomselect top "$atsel_in" frame $frame_i]
    set lipid_com_z [lindex [measure center $lipidsel weight mass] 2]

    if {$lipid_com_z < $refsel_com_z} {
        $lipidsel set user2 -1
        $lipidsel delete
        return -1
    } else {
        $lipidsel set user2 1
        $lipidsel delete
        return 1
    }
}

;# originally by Grace Brannigan, named local_midplane2
;# Selects all PO4/GL1/GL2/AM1/AM2 within ~1.4 nm of lipid's COM x,y coordinate.
;# Interprets COM of PO4/GL1/GL2/AM1/AM2 selection z component to be the "local midplane."
;# Compares lipid's COM z component to local midplane and sorts accordingly.
proc leaflet_sorter_3 {atsel_in frame_i} {
    set lipidsel [atomselect top $atsel_in frame $frame_i]
    set lipid_com [measure center $lipidsel weight mass]
    set lipid_x [lindex $lipid_com 0]
    set lipid_y [lindex $lipid_com 1]
    set lipid_z [lindex $lipid_com 2]
    
    set local_surfaces [atomselect top "name PO4 GL1 GL2 AM1 AM2 and pbwithin 200 of $atsel_in" frame $frame_i]
    set local_midplane [lindex [measure center $local_surfaces weight mass] 2]
    $local_surfaces delete

    if {$lipid_z < $local_midplane} {
        $lipidsel set user2 -1
        $lipidsel delete
        return -1
    } else {
        $lipidsel set user2 1
        $lipidsel delete
        return 1
    }
}

;# Determines if the lipid is in the outer or inner leaflet and sets the user2 value accordingly
;# Algorithm is selected by user specifying one of the following options: 
;# 0: determines leaflet based on relative height of specified head and tail beads
;# 1: originally by Liam Sharp; procedure that was used in JCP 2021 for nAChR; similar to leaflet_sorter_0 but autoselects head and tail beads; more appropriate for situations with many species
;# 2: originally by Jahmal Ennis; determines whether the auto-determined headbead is above or below the center of mass of some reference selection; more appropriate for rigid lipids like cholesterol that frequently invert or lie at parallel to the membrane
;# 3: originally by Grace Brannigan and called local_midplane2; selects all PO4/GL1/GL2 beads within a circular region around the lipid, measures the COM, assumes that COM to be the midplane, and sorts lipids based on whether their COM is above or below the midplane.
;# 4: assigns all molecules to the outer leaflet
proc leaflet_detector {atsel_in head tail frame_i leaflet_sorting_algorithm} {
    global params
    if {$leaflet_sorting_algorithm == 0} {
        leaflet_sorter_0 $atsel_in $head $tail $frame_i
    } elseif { $leaflet_sorting_algorithm == 1 } {
        leaflet_sorter_1 $atsel_in $frame_i
    } elseif { $leaflet_sorting_algorithm == 2 } {
        leaflet_sorter_2 $atsel_in $params(leaflet_sorter_2_reference_sel) $frame_i
    } elseif { $leaflet_sorting_algorithm == 3 } {
        leaflet_sorter_3 $atsel_in $frame_i
    } elseif { $leaflet_sorting_algorithm == 4} {
        set sel [atomselect top $atsel_in frame $frame_i]
        $sel set user2 1
        $sel delete
    } else { 
        #default
        leaflet_sorter_1 $atsel_in $frame_i
    }
}


;# Calculates the total number of lipids and beads of the given selection in each leaflet 
;# Assigns the leaflet to user2 
;# Returns the following list : [["lower" lower_leaflet_beads lower_leaflet_lipids] ["upper" upper_leaflet_beads upper_leaflet_lipids]] 
proc frame_leaflet_assignment {atseltext headname tailname frame_i frame_f {restrict_to_Rmax 0}} {
    global params
    if {$restrict_to_Rmax == 1} {
        set outer_r2 [expr $params(Rmax)**2]
        set sel [ atomselect top "($atseltext) and same resid as (x*x + y*y < $outer_r2)"  frame $frame_i]
    } elseif {$restrict_to_Rmax == 0} {
        set sel [ atomselect top "$atseltext"  frame $frame_i]
    } else {
        error "restrict_leaflet_sorter_to_Rmax must be 1 or 0"
    }
    set sel_num [llength [lsort -unique [$sel get resid] ] ]
    set sel_resid_list [lsort -unique [$sel get resid] ]
    set totals {}
    if {$sel_num < 1} {
        set totals [list "lower 0 0" "upper 0 0"] 
    } else {
        #assign leaflets from $frame_i to user2 field of each bead for this selection
        foreach sel_resid $sel_resid_list {
            set selstring "(${atseltext}) and (resid $sel_resid)"
            leaflet_detector $selstring $headname $tailname $frame_i $params(leaflet_sorting_algorithm)
        }
        #copy leaflet values from $frame_i to all frames between $frame_i and $frame_f
        #need a new atomselection that is just based on index rather than shell
        set leaflet_list [$sel get user2] 
        set indices [$sel get index]
        set interim_sel [atomselect top "index $indices"]
        for {set interim_frame [expr $frame_i + 1]} {$interim_frame < [expr $frame_f]} {incr interim_frame} {
            $interim_sel frame $interim_frame
            $interim_sel update
            $interim_sel set user2 $leaflet_list
        }
        $interim_sel delete
        #count the number of lipids and the number of beads in each leaflet
        foreach leaf [list  "(user2<0)" "(user2>0)"] txtstr [list "lower" "upper"] {
            set leaf_sel [ atomselect top "(${atseltext}) and $leaf"  frame $frame_i]
            set num_beads [$leaf_sel num]
            set num_lipids [llength [lsort -unique [$leaf_sel get resid] ]]
            lappend totals [list $txtstr $num_beads $num_lipids]
            $leaf_sel delete
        }
    }
    $sel delete
    return $totals
}

;# Calculates the total number of lipids and beads of the given selection in each leaflet 
;# Returns the following list : [["lower" lower_leaflet_beads lower_leaflet_lipids] ["upper" upper_leaflet_beads upper_leaflet_lipids]] 
proc trajectory_leaflet_assignment {atseltext headname tailname} { 
    global params
    set num_reassignments 0
    if {[lsearch -exact "0 1 2 3 4" $params(leaflet_sorting_algorithm)] == -1} {
        puts "Option $params(leaflet_sorting_algorithm) not recognized as a leaflet sorting option. Defaulting to option 1."
    } elseif {$params(leaflet_sorting_algorithm) == 2} {
        if {$params(leaflet_sorter_2_reference_sel) eq "none"} {
                puts "No reference selection provided for leaflet sorter 2."
                puts "Defaulting to z=0 as the reference height to sort by."
        }
    }
    for {set update_frame $params(start_frame)} {$update_frame <= [expr $params(end_frame) - $params(leaflet_reassign_interval)]} {incr update_frame $params(leaflet_reassign_interval)} {
        frame_leaflet_assignment $atseltext $headname $tailname $update_frame [expr $update_frame + $params(leaflet_reassign_interval)] $params(restrict_leaflet_sorter_to_Rmax)
        incr num_reassignments
    }
    if {[test_if_evenly_divisible $params(end_frame) $params(leaflet_reassign_interval)] != 1} {
        # Run one extra iteration to finish final leftover frames at end of trajectory
        frame_leaflet_assignment $atseltext $headname $tailname $update_frame $params(end_frame) $params(restrict_leaflet_sorter_to_Rmax)
        incr num_reassignments
    }
    puts "Checked for leaflet reassignments $num_reassignments times."
}
    
;#Reinitializes the user2 value for selected beads in selected frames 
proc clean_leaflet_assignments {atseltext} {
    global params
    set sel [ atomselect top "$atseltext"]
    set selnum [$sel num]

    for {set update_frame $params(start_frame)} {$update_frame < ${ _frame}} {incr update_frame} {
        $sel frame $update_frame
        $sel set user2 [lrepeat $selnum 0.0]
        puts "Cleaning $selnum beads of leaflet assignments in frame $update_frame"
    }
    $sel delete
}
    
#test to see if two floats are evenly divisible. Return 1 if evenly divisible.
#Return 0 if not evenly divisible.
proc test_if_evenly_divisible {dividend divisor} {
    set TOLERANCE [expr 10.0**-12]
    set float_quotient [expr $dividend / double($divisor)]
    set int_quotient [expr int($float_quotient)]
    set diff [expr $float_quotient - $int_quotient]
    if {$diff <= $TOLERANCE} {
        return 1
    } else {
        return 0
    }
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


;#The inner-most loop of the histogramming algorithm: a loop over all lipid atoms (or beads) occupying one shell in one frame. Each atom is assigned an angular bin and totals are updated.  
proc loop_over_atoms {shell atseltext frm} {
    global params
    set indexs [$shell get index]
    set theta_high_out [list]
    set theta_low_out [list]
    set leaflet 0
    foreach indx $indexs {
        #loop over atoms (or beads if CG) in the shell
        set atsel "($atseltext) and index $indx"
        set thislipid [atomselect top $atsel frame $frm]
        set x [$thislipid get x]
        set y [$thislipid get y]
        set leaflet [$thislipid get user2]
        set theta [get_theta $x $y]
        set theta_bin [expr int($theta/$params(dtheta))]
        if {$leaflet > 0} {
            lappend theta_high_out $theta_bin
        } elseif {$leaflet < 0} {
            lappend theta_low_out $theta_bin
        } else {
            puts "WARNING: lipid atom $indx did not get assigned a leaflet for frame $frm"
        }
        $thislipid set user [expr $theta_bin]
        $thislipid delete
    }
    
    return [list $theta_low_out $theta_high_out] ;#lower before upper is the convention
}

;#The middle nested loop of the histogramming algorithm: a loop over all frames for a given radial shell. The atoms/beads occupying the shell are calculated using atomselect within and updated in each frame, without creating or destroying a new atom selection. 
proc loop_over_frames {shell atseltext start_frame end_frame ri rf flower fupper r_index} {
    global params
    set theta_bin_high [lrepeat [expr $params(Ntheta)+1] 0]
    set theta_bin_low [lrepeat [expr $params(Ntheta)+1] 0]
    for {set frm $params(start_frame)} {$frm < ${end_frame}} {incr frm $params(dt)} {
        $shell frame $frm
        $shell update 
        $shell set user3 $r_index
        set singleFrame_counts [loop_over_atoms $shell $atseltext $frm]
        set singleFrame_upper [lindex $singleFrame_counts 1] 
        set singleFrame_lower [lindex $singleFrame_counts 0]
        set theta_bins [theta_histogram $singleFrame_lower $singleFrame_upper]
        if { [llength $theta_bin_high] != [llength [lindex $theta_bins 0]] } {
            error "theta_bin_high/low and theta_bins do not have the same length."
        }
        set theta_bin_high [vecadd $theta_bin_high [lindex $theta_bins 1] ]
        #puts [lindex $theta_bins 1]
        set theta_bin_low [vecadd $theta_bin_low [lindex $theta_bins 0]]
        #puts $theta_bin_low
        output_bins $fupper $ri $rf [lindex $theta_bins 1] 
        output_bins $flower $ri $rf [lindex $theta_bins 0]   
    }
    return [list  ${theta_bin_low} ${theta_bin_high}]
}


;#The outer loop of the 3 nested histogramming loops. 
;#Unintuitively, the outermost loop is over radial shells, then the middle loop is over frames, and the inner most loop is over atoms/beads in the shell. 
;#This odd construction improves efficiency: the radial atomselections can be created using "atomselect within", and then updated for each new frame in the middle loop, without a new selection being created or destroyed. There is no equivalent option for an angular "within" so angular histogramming occurs more traditionally via a loop over atoms/beads.  

proc loop_over_shells {atseltext low_f upp_f low_f_avg upp_f_avg} {
    global params
    set delta_frame [expr ($params(end_frame) - $params(start_frame)) / $params(dt)]
    set radial_bin_index 0
    for {set ri $params(Rmin)} { $ri<$params(Rmax)} { set ri [expr $ri + $params(dr)]} {
        #loop over shells
        puts "Now on shell {$ri [expr ${ri}+$params(dr)]}"
        set rf [expr $ri + $params(dr)]
        set rf2 [expr $rf*$rf]
        set ri2 [expr $ri*$ri]
        set shell [atomselect top "($atseltext) and ((x*x + y*y < $rf2) and  (x*x + y*y > $ri2))"]
        set theta_bin [loop_over_frames $shell $atseltext $params(start_frame) $params(end_frame) $ri $rf $low_f $upp_f $radial_bin_index]
        set theta_bin_high [lindex $theta_bin 1]
        set theta_bin_low [lindex $theta_bin 0]
        $shell delete	
        set time_avg_upper [vecscale $theta_bin_high [expr 1.0 / (1.0 * $delta_frame)]]
        set time_avg_lower [vecscale $theta_bin_low [expr 1.0 / (1.0 * $delta_frame)]]
        output_bins $upp_f_avg $ri $rf "$time_avg_upper" 
        output_bins $low_f_avg $ri $rf "$time_avg_lower" 
        incr radial_bin_index
    }
}


proc set_parameters { config_file_script } {
    global params
    array unset params
    global params
    array set params {
        leaflet_sorting_algorithm 1
        leaflet_sorter_2_reference_sel "none"
        center_and_align 0
        use_qwrap 0
        utils ./helpers 
        dt 1
        leaflet_reassign_interval 1
        start_frame 0  
        backbone_selstr "name BB" 
        protein_selstr "name BB SC1 to SC4"
        atomsels {"resname POPG"}
        headnames {"name PO4"}
        tailnames {"name C4A C4B"}
        chainlist {A B C D E}
        helixlist {1 2 3 4}
        helix_assignment_script "assign_helices_ELIC_general.tcl" 
        midplane_selstr "occupancy 1 to 4"        
        Rmax 20.
        Rmin 0.
        dr 5 
        Ntheta 50
        restrict_leaflet_sorter_to_Rmax 0
        filename_stems {"POPG"}
    }
    set nframes [molinfo top get numframes] 
    array set params [list end_frame $nframes]

    set param_name_list [lsort -dictionary [array names params]]
    source $config_file_script
    set unread_params {}
    puts "------------Configuration Parameters------------"
    foreach param_name $param_name_list {
        if {[info exists $param_name]} {
            set param_value [set $param_name] 
            array set params [list $param_name $param_value]
            puts "$param_name: $params($param_name)"
        } else {
            lappend unread_params $param_name
        }
    }
    foreach param_name $unread_params {
        puts "Warning: No value of parameter $param_name was read from the config file. Using the default value of $params($param_name)."         
    }
    puts "------------------------------------------------"
    array set params [list dtheta [expr 360.0/(1.0*($params(Ntheta)))]]
    return params
}





### polarDensity Function ###


;#The main function that initializes, constructs the densities for each lipid selection, and outputs to file. 

proc polarDensityBin { config_file_script } { 
    ;#read parameters
    global params
    set_parameters $config_file_script
    source $params(utils)/BinTools.tcl

    ;# check to make sure Rmax is evenly divisible by dr.
    if {[test_if_evenly_divisible $params(Rmax) $params(dr)] != 1} {
        error "Rmax must be evenly divisible by dr."
    }
    if {$params(use_qwrap) == 1} {load $params(utils)/qwrap.so}
    set backbone_selstr $params(backbone_selstr) ;#only necessary for backwards compatibility 
    set protein_selstr $params(protein_selstr) ;#only necessary for backwards compatibility 
    source $params(helix_assignment_script)
    foreach atseltext $params(atomsels) stem $params(filename_stems) headname $params(headnames) tailname $params(tailnames) {
        ;# make sure the atomselection exists
        set sel [atomselect top "$atseltext"]
        set sel_num [$sel num]
        $sel delete
        if {$sel_num == 0} {
            error "No lipid of matching atomselection text $atseltext"
        }

        if {$params(center_and_align) == 1} {
            ;# wraps and centers system at origin
            center_and_wrap_system "occupancy $params(helixlist) and $params(backbone_selstr)"
            ;# aligns protein
            Align "occupancy $params(helixlist) and $params(backbone_selstr)"
        }
        ;# outputs protein positions
        output_inclusion_centers
        ;# initialize some constants
        set area [get_avg_area top]
        set nframes [molinfo top get numframes]
        if { $params(start_frame) > $nframes } {
            puts "Error: specified start frame $params(start_frame) is greater than number of frames $nframes" 
            set start_frame $nframes
        }
        if { $params(end_frame) > $nframes } {
            puts "Warning: specified end frame $params(end_frame) is greater than number of frames; setting end frame to $nframes" 
            set end_frame $nframes
        }

        puts "Atomselection:\t$atseltext"
        set low_f [open "${stem}.low.dat" w]
        set upp_f [open "${stem}.upp.dat" w]
        set low_f_avg [open "${stem}.low.avg.dat" w]
        set upp_f_avg [open "${stem}.upp.avg.dat" w]
        set totals [frame_leaflet_assignment $atseltext $headname $tailname $params(end_frame) $params(end_frame)]   
        
        foreach lu [list $low_f $upp_f] avgfile [list $low_f_avg $upp_f_avg] leaf_total $totals {
            set leaflet_str [lindex $leaf_total 0]
            set expected_beads [lindex $leaf_total 1]
            set expected_lipids [lindex $leaf_total 2]
            set expected_bead_density [expr 1.0 * $expected_beads/$area]
                puts "#Selection name $stem in $leaflet_str leaflet: ${expected_lipids} molecules, Num beads : ${expected_beads} beads,  Average Area : [format {%0.0f} $area] A^2, Expected Bead Density : [format {%0.5f} [expr $expected_bead_density]]/A^2, dr*dtheta : [format {%0.5f} [expr $params(dr)*[DtoR $params(dtheta)]]] "
                puts $lu "#Selection name $stem in $leaflet_str leaflet: ${expected_lipids} molecules, Num beads : ${expected_beads} beads,  Average Area : [format {%0.0f} $area] A^2, Expected Bead Density : [format {%0.5f} [expr $expected_bead_density]]/A^2, dr*dtheta : [format {%0.5f} [expr $params(dr)*[DtoR $params(dtheta)]]] "
                puts $avgfile "#Selection name $stem in $leaflet_str leaflet: ${expected_lipids} molecules, Num beads : ${expected_beads} beads,  Average Area : [format {%0.0f} $area] A^2, Expected Bead Density : [format {%0.5f} [expr $expected_bead_density]]/A^2, dr*dtheta : [format {%0.5f} [expr $params(dr)*[DtoR $params(dtheta)]]] "
        }
        puts "Processing frames, starting at frame $params(start_frame) and ending at frame $params(end_frame)."
        trajectory_leaflet_assignment $atseltext $headname $tailname   
        ;#the core calculation 
        loop_over_shells $atseltext $low_f $upp_f $low_f_avg $upp_f_avg  
        
        close $low_f
        close $upp_f
        close $low_f_avg
        close $upp_f_avg
    }
}    