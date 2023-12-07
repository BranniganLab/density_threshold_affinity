
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
    set pi 3.14159265358979323846
    return [expr $r*180.0/$pi]
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


# z_mid
#
# Finds the average mid-plane
# Arguments:
#   int: first frame
#   nframes: The number of frames over which to average
# Outputs:
#   float: The average z value of all the beads
#

proc z_mid {init_frm nframes midplane_selstr} {
    set z_list {}
    for {set frm ${init_frm}} {${frm} < ${nframes}} {incr frm} {
        set mid [atomselect top $midplane_selstr frame $frm]
        lappend z_list [lindex [measure center $mid weight mass] 2]
        $mid delete
    }
    return [expr 1.0*[vecsum $z_list]/([llength $z_list]) ]
}



;#Outputs xy position of helix centers to file for each leaflet; center is calculated using the part of the helix in the given leaflet
proc output_helix_centers {chain_names helix_occupancy_list  backbone_selstr midplane_selstr {a ""} } {
    ;# list for the chain names
    ;# finds the center of the membranes
    set zed [z_mid 0 20 $midplane_selstr]
    ;# calculates the center of mass for subunit alpha helices in both leaflets
    puts "Writing coordinates for [llength $chain_names] chains and [llength $helix_occupancy_list] helices per chain"
    foreach eq {"<" ">"} eqtxt {"lwr" "upr"} {
        set fout [open "./Protein${a}_coords_${eqtxt}.dat" w]
        puts $fout  "# chain A ooc1r occ1the occ2r occ2the... "
        foreach chnm $chain_names {
            foreach occ $helix_occupancy_list {
                set sel [atomselect top "(chain ${chnm}) and (occupancy $occ and $backbone_selstr) and (z ${eq} $zed)" frame 0]
                set com [measure center $sel weight mass]
                $sel delete
                set x [lindex $com 0]
                set y [lindex $com 1]
                set r [expr sqrt($x*$x+$y*$y)]
                set theta [get_theta $x $y]
                #puts "chain ${chnm} and occupancy $occ $r $theta"
                puts -nonewline $fout "$r $theta "
            }
            puts $fout ""
        }
        close $fout
    }
}

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


proc center_system {inpt} {
    global USE_QWRAP
    puts "${inpt}"
    # confirms your box is either square or paraelleogram-ish
    # will perform qwrap or pbc wrap depending
    
    set pbc_angles [molinfo top get {alpha beta gamma}]
    
    set sel [atomselect top "$inpt" frame 0]
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
        
        if {($USE_QWRAP==0) || (([lindex $pbc_angles 0]!=90.0) && ([lindex $pbc_angles 0]!=90.0) && ([lindex $pbc_angles 0]!=90.0))} {
            puts "qwrap may not be optimal for your system...\n"
            puts "Running pbc wrap. To verify proper centering"
            puts "pbc wrap will be run multiple times" ; after 100
            foreach i {0 1 2 3} {
                pbc wrap -centersel "$inpt" -all
            }
        } else {
            qwrap centersel "$inpt" ;#center entire system at ~0,0,0
        }
        set com [measure center $sel weight mass]
        incr counter_i
    }
    $sel delete
}

proc leaflet_sorter_0 {atsel_in head tail frame_i} {
    ;#originally by Grace Brannigan
    #puts "Sorting into leaflets using leaflet_sorter_0"
    set sel_resid [atomselect top "$atsel_in" frame $frame_i]
    set sel_head [atomselect top "$atsel_in and name $head" frame $frame_i]
    set sel_tail [atomselect top "$atsel_in and name $tail" frame $frame_i]
    
    set head_Z [${sel_head} get z] 
    set tail_Z [${sel_tail} get z] 
    
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


proc leaflet_sorter_1 {atsel_in frame_i} {
    ;#originally by Liam Sharp; procedure that was used in JCP 2021 for nAChR
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


proc leaflet_sorter_2 {atsel_in frame_i} {
    ;#originally by Jahmal Ennis, designed for cholesterol 
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


;# Determines if the lipid is in the outer or iner leaflet and sets the user value accordingly
;# Returns +1 if the lipid is in the upper leaflet and -1 if it is in the lower leaflet 
proc leaflet_detector {atsel_in head tail frame_i leaflet_algorithm} {
    if {$leaflet_algorithm == 0} {
        leaflet_sorter_0 $atsel_in $head $tail $frame_i
    } elseif { $leaflet_algorithm == 1 } {
        leaflet_sorter_1 $atsel_in $frame_i
    } elseif { $leaflet_algorithm == 2 } {
        leaflet_sorter_2 $atsel_in $frame_i
    } else { 
        puts "Option $LEAFLET_SORTING_ALGORITHM not recognized as a leaflet sorting option.  Defaulting to option 1." 
    }
}

;# Unnecessary; duplicated by frame_leaflet_assignment with frame_f = frame_i
;# Calculates the total number of lipids and beads of the given species in each leaflet 
;# Returns the following list : [["lower" lower_leaflet_beads lower_leaflet_lipids] ["upper" upper_leaflet_beads upper_leaflet_lipids]]  
proc get_leaflet_totals {species headname tailname lipidbeads_selstr frame_i leaflet_algorithm} {
    set sel [ atomselect top "(($species)) and $lipidbeads_selstr"  frame $frame_i]
    set sel_num [llength [lsort -unique [$sel get resid] ] ]
    set sel_resid_list [lsort -unique [$sel get resid] ]
    set totals {}
    $sel delete
    if {$sel_num < 1} {
        set totals [[list 0 0] [list 0 0]] 
    } else {
        #assign leaflets to user2 field of each bead for this species
        foreach sel_resid $sel_resid_list {
            set selstring "${species} and (resid $sel_resid) and $lipidbeads_selstr"
            set leaflet [leaflet_detector $selstring $headname $tailname $frame_i $leaflet_algorithm]
        }   
        #count the number of lipids and the number of beads in each leaflet
        foreach leaf [list  "(user2<0)" "(user2>0)"] txtstr [list "lower" "upper"] {
            set sel [ atomselect top "(${species} and $leaf)"  frame $frame_i]
            set num_beads [$sel num]
            set num_lipids [llength [lsort -unique [$sel get resid] ]]
            lappend totals [list $txtstr $num_beads $num_lipids]
            $sel delete
        }
    }
    return $totals
}

;# Calculates the total number of lipids and beads of the given species in each leaflet 
;# Assigns the leaflet to user2 
;# Returns the following list : [["lower" lower_leaflet_beads lower_leaflet_lipids] ["upper" upper_leaflet_beads upper_leaflet_lipids]] 
proc frame_leaflet_assignment {species headname tailname lipidbeads_selstr frame_i frame_f leaflet_algorithm} {
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
            set leaflet [leaflet_detector $selstring $headname $tailname $frame_i $leaflet_algorithm]
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
;# Returns the following list : [[lower_leaflet_beads lower_leaflet_lipids] [upper_leaflet_beads upper_leaflet_lipids]] 
proc trajectory_leaflet_assignment {species headname tailname lipidbeads_selstr start end skip leaflet_algorithm} {
    set num_reassignments 0
    for {set update_frame $start} {$update_frame < ${end}} {incr update_frame $skip} {
        frame_leaflet_assignment $species $headname $tailname $lipidbeads_selstr $update_frame [expr $update_frame + $skip] $leaflet_algorithm
        incr num_reassignments
    }
    puts "Checked for leaflet reassignments $num_reassignments times."
}
    
;# Calculates the total number of lipids and beads of the given species in each leaflet 
;# Returns the following list : [[lower_leaflet_beads lower_leaflet_lipids] [upper_leaflet_beads upper_leaflet_lipids]] 
proc clean_leaflet_assignments {species lipidbeads_selstr start end} {
    set sel [ atomselect top "$species and $lipidbeads_selstr"]
    set selnum [$sel num]

    for {set update_frame $start} {$update_frame < ${end}} {incr update_frame} {
        $sel frame $update_frame
        $sel set user2 [lrepeat $selnum 0.0]
        puts "Cleaning $selnum beads of leaflet assignments in frame $update_frame"
    }
    $sel delete
}
    
    
proc output_bins {fl  ri rf dtheta bins} {
    puts -nonewline $fl "[format {%0.2f} $ri]  [format {%0.2f} $rf] [format {%0.2f} $dtheta]  " 
    puts $fl "$bins" 
}

proc theta_histogram {singleFrame_lower singleFrame_upper  Ntheta } {
    
    set theta_bin_out [list]
    
    foreach ud [list $singleFrame_lower $singleFrame_upper ] {
        #cleanup and output 
        set theta_bin_counts [lcount $ud]
        #Shell_Test $shel_count $theta_bin_counts
        set theta_bins {}
        # make this into the new lcount? better Idea TEST lcount
        for {set ti 0} { $ti<=$Ntheta} {incr ti 1} {
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


#bins a single shell over many frames
proc bin_over_frames {shell species headname tailname lipidbeads_selstr dtheta start end Ntheta dt ri rf  flower fupper leaflet_algorithm} {
    set theta_bin_high [lrepeat [expr $Ntheta+1] 0]
    set theta_bin_low [lrepeat [expr $Ntheta+1] 0]
    for {set frm $start} {$frm < ${end}} {incr frm $dt} {
        $shell frame $frm
        $shell update 
        set singleFrame_counts [bin_frame $shell $species $headname $tailname $lipidbeads_selstr $dtheta $frm $leaflet_algorithm ]
        set singleFrame_upper [lindex $singleFrame_counts 1] 
        set singleFrame_lower [lindex $singleFrame_counts 0]
        set theta_bins [theta_histogram $singleFrame_lower $singleFrame_upper  $Ntheta]
        
        # should be fixed, do not change [lrepeat [expr $Ntheta+1] to [lrepeat [expr $Ntheta] 
        if { [llength $theta_bin_high] != [llength [lindex $theta_bins 0]] } {
            error "theta_bin_high/low and theta_bins do not have the same length."
        }
        set theta_bin_high [vecadd $theta_bin_high [lindex $theta_bins 1] ]
        #puts [lindex $theta_bins 1]
        set theta_bin_low [vecadd $theta_bin_low [lindex $theta_bins 0]]
        #puts $theta_bin_low
        output_bins $fupper $ri $rf $dtheta [lindex $theta_bins 1] 
        ;#open fupper before the loop starts and close afterwards
        output_bins $flower $ri $rf $dtheta [lindex $theta_bins 0] 
        ;#same thing     
        
    }
    return [list  ${theta_bin_low} ${theta_bin_high}]
}


proc bin_frame {shell species headname tailname lipidbeads_selstr dtheta frm leaflet_algorithm} {
    set indexs [$shell get index]
    set resids [$shell get resid]
    set nShell [$shell num]
    set theta_high_out [list]
    set theta_low_out [list]
    set resd_old 0
    set leaflet 0
    foreach indx $indexs resd $resids {
        #loop over lipids in the shell
        set a "($species and index $indx)"
        set b "(resid $resd and $species and $lipidbeads_selstr)" 
        set thislipid [atomselect top $a frame $frm]
#       if {[string length ${species}] == 2} {
#           if {([$thislipid get name] == "PO4") || ([$thislipid get name] == "P") } { ;#GB has no idea what this does. 
#               continue
#           }
#       }
        set x [$thislipid get x]
        set y [$thislipid get y]
        set leaflet [$thislipid get user2] ;
        set theta [get_theta $x $y]
        set ti [expr int($theta/$dtheta)] 
        if {$leaflet > 0} {
            lappend theta_high_out $ti
        } elseif {$leaflet <0} {
            lappend theta_low_out $ti
        } else {
            puts "WARNING: lipid $resd did not get assigned a leaflet for frame $frm"
        }
        $thislipid set user [expr $ti+1]
        $thislipid delete
    }
    
    return [list $theta_low_out $theta_high_out] ;#lower before upper is the convention
}





### polarDensity Function ###



proc polarDensityBin { config_file_script } { 
    set start_frame 0 ; #default value before potential change in $config_file_script
    set nframes [molinfo top get numframes]
    set end_frame $nframes ;#default value before potential change in $config_file_script
    source $config_file_script
    source $UTILS/BinTools.tcl
    if {$USE_QWRAP == 1} {load ${UTILS}/qwrap.so}
    source ${helix_assignment_script}
    
    foreach species $lipids lipidbeads_selstr $lipidbeads_selstrs acylchain_selstr $acylchain_selstrs headname $headnames tailname $tailnames {
        set outfile "$species"
        set sel [atomselect top "resname $species"]
        set sel_num [$sel num]
        
        if {$sel_num == 0} {
            error "No lipid of species $species"
        }
        ;# Center's system (weak hack)
        if {$CENTER_AND_ALIGN == 1} {
            Center_System "occupancy $helixlist and $backbone_selstr"
            Center_System "occupancy $helixlist and $backbone_selstr"
            Center_System "occupancy $helixlist and $backbone_selstr"
            ;# aligns protein
            Align "occupancy $helixlist and $backbone_selstr"
        }
        ;# outputs protein positions
        output_helix_centers $chainlist $helixlist $backbone_selstr $midplane_selstr
        ;# initialize some constants
        set area [get_avg_area top]

        if { $start_frame > $nframes } {
            puts "Error: specified start frame $start_frame is greater than number of frames $nframes" 
            set end $nframes
        }
        if { $end_frame > $nframes } {
            puts "Warning: specified end frame $end_frame is greater than number of frames; setting end frame to $nframes" 
            set end $nframes
        }
        $sel delete
        puts "Acyl Chain:\t$species"
        set low_f [open "${outfile}.low.dat" w]
        set upp_f [open "${outfile}.upp.dat" w]
        set low_f_avg [open "${outfile}.low.avg.dat" w]
        set upp_f_avg [open "${outfile}.upp.avg.dat" w]
        set dtheta [expr 360.0/(1.0*($Ntheta))]
        set totals [frame_leaflet_assignment "resname $species" $headname $tailname $lipidbeads_selstr $start_frame $start_frame $LEAFLET_SORTING_ALGORITHM ]
        
        foreach lu [list $low_f $upp_f] avgfile [list $low_f_avg $upp_f_avg] leaf_total $totals {
            set leaflet_str [lindex $leaf_total 0]
            set expected_beads [lindex $leaf_total 1]
            set expected_lipids [lindex $leaf_total 2]
            set expected_bead_density [expr 1.0 * $expected_beads/$area]
                puts "#Lipid species $species in $leaflet_str leaflet: ${expected_lipids} molecules, Num beads : ${expected_beads} beads,  Average Area : [format {%0.0f} $area] A^2, Expected Bead Density : [format {%0.5f} [expr $expected_bead_density]]/A^2, Average Chain : [avg_acyl_chain_len "resname $species" $acylchain_selstr] beads, dr*dtheta : [format {%0.5f} [expr $dr*[DtoR $dtheta]]] "
                puts $lu "#Lipid species $species in $leaflet_str leaflet: ${expected_lipids} molecules, Num beads : ${expected_beads} beads,  Average Area : [format {%0.0f} $area] A^2, Expected Bead Density : [format {%0.5f} [expr $expected_bead_density]]/A^2, Average Chain : [avg_acyl_chain_len "resname $species" $acylchain_selstr] beads, dr*dtheta : [format {%0.5f} [expr $dr*[DtoR $dtheta]]] "
                puts $avgfile "#Lipid species $species in $leaflet_str leaflet: ${expected_lipids} molecules, Num beads : ${expected_beads} beads,  Average Area : [format {%0.0f} $area] A^2, Expected Bead Density : [format {%0.5f} [expr $expected_bead_density]]/A^2, Average Chain : [avg_acyl_chain_len "resname $species" $acylchain_selstr] beads, dr*dtheta : [format {%0.5f} [expr $dr*[DtoR $dtheta]]] "
        }
        
        trajectory_leaflet_assignment "resname $species" $headname $tailname $lipidbeads_selstr $start_frame $end_frame $leaflet_reassign_t $LEAFLET_SORTING_ALGORITHM
        set delta_frame [expr ($end_frame - $start_frame) / $dt]
        for {set ri $Rmin} { $ri<=${Rmax}} { set ri [expr $ri + $dr]} {
            #loop over shells
            puts "Now on shell {$ri [expr ${ri}+${dr}]}"
            set rf [expr $ri + $dr]
            set rf2 [expr $rf*$rf]
            set ri2 [expr $ri*$ri]
            set shell [atomselect top "(resname $species) and ((x*x + y*y < $rf2) and  (x*x + y*y > $ri2)) and $lipidbeads_selstr"]
            #puts [$shell num]		
            set theta_bin [bin_over_frames $shell "resname $species" $headname $tailname $lipidbeads_selstr $dtheta $start $end $Ntheta $dt $ri $rf $low_f $upp_f $LEAFLET_SORTING_ALGORITHM]
            set theta_bin_high [lindex $theta_bin 1]
            set theta_bin_low [lindex $theta_bin 0]
            $shell delete	
            set time_avg_upper [vecscale $theta_bin_high [expr 1.0 / (1.0 * $delta_frame)]]
            set time_avg_lower [vecscale $theta_bin_low [expr 1.0 / (1.0 * $delta_frame)]]
            output_bins $upp_f_avg $ri $rf $dtheta "$time_avg_upper" 
            output_bins $low_f_avg $ri $rf $dtheta "$time_avg_lower" 
        }
        close $low_f
        close $upp_f
        close $low_f_avg
        close $upp_f_avg
    }
}    