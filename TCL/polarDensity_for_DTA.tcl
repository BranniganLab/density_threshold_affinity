
package require pbctools

;#Outputs xy position of helix centers to file for each leaflet; center is calculated using the part of the helix in the given leaflet
proc output_helix_centers {{a ""} } {
    global params
    ;# list for the chain names
    ;# finds the center of the membranes
    set zed [binningTools::z_mid 0 20]
    puts "Liam Says : $zed"
    ;# calculates the center of mass for subunit alpha helices in both leaflets
    puts "Writing coordinates for [llength $params(chainlist)] chains and [llength $params(helixlist)] helices per chain"
    foreach eq {"<" ">"} eqtxt {"lwr" "upr"} {
        set fout [open "./Protein${a}_coords_${eqtxt}.dat" w]
        puts $fout  "# chain A ooc1r occ1the occ2r occ2the... "
        foreach chnm $params(chainlist) {
            foreach occ $params(helixlist) {
                set sel [atomselect top "(chain ${chnm}) and (occupancy $occ and $params(backbone_selstr)) and (z ${eq} $zed)" frame 0]
                set com [measure center $sel weight mass]
                $sel delete
                set x [lindex $com 0]
                set y [lindex $com 1]
                set r [expr sqrt($x*$x+$y*$y)]
                set theta [mathTools::get_theta $x $y]
                #puts "chain ${chnm} and occupancy $occ $r $theta"
                puts -nonewline $fout "$r $theta "
            }
            puts $fout ""
        }
        close $fout
    }
}

;#The inner-most loop of the histogramming algorithm: a loop over all lipids occupying one shell in one frame. Each lipid is assigned an angular bin and totals are updated.  
proc loop_over_lipids {shell species headname tailname lipidbeads_selstr frm} {
    global params
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
        set x [$thislipid get x]
        set y [$thislipid get y]
        set leaflet [$thislipid get user2] ;
        set theta [::mathTools::get_theta $x $y]
        set ti [expr int($theta/$params(dtheta))] 
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

;#The middle nested loop of the histogramming algorithm: a loop over all frames for a given radial shell. The lipids occupying the shell are calculated using atomselect within and updated in each frame, without creating or destroying a new atom selection. 
proc loop_over_frames {shell species headname tailname lipidbeads_selstr start_frame end_frame ri rf flower fupper} {
    global params
    set theta_bin_high [lrepeat [expr $params(Ntheta)+1] 0]
    set theta_bin_low [lrepeat [expr $params(Ntheta)+1] 0]
    for {set frm $params(start_frame)} {$frm < ${end_frame}} {incr frm $params(dt)} {
        $shell frame $frm
        $shell update 
        set singleFrame_counts [loop_over_lipids $shell $species $headname $tailname $lipidbeads_selstr $frm]
        set singleFrame_upper [lindex $singleFrame_counts 1] 
        set singleFrame_lower [lindex $singleFrame_counts 0]
        set theta_bins [binningTools::theta_histogram $singleFrame_lower $singleFrame_upper]
        if { [llength $theta_bin_high] != [llength [lindex $theta_bins 0]] } {
            error "theta_bin_high/low and theta_bins do not have the same length."
        }
        set theta_bin_high [vecadd $theta_bin_high [lindex $theta_bins 1] ]
        #puts [lindex $theta_bins 1]
        set theta_bin_low [vecadd $theta_bin_low [lindex $theta_bins 0]]
        #puts $theta_bin_low
        binningTools::output_bins $fupper $ri $rf [lindex $theta_bins 1] 
        binningTools::output_bins $flower $ri $rf [lindex $theta_bins 0]   
    }
    return [list  ${theta_bin_low} ${theta_bin_high}]
}


;#The outer loop of the 3 nested histogramming loops. 
;#Unintuitively, the outermost loop is over radial shells, then the middle loop is over frames, and the inner most loop is over lipids in the shell. 
;#This odd construction improves efficiency: the radial atomselections can be created using "atomselect within", and then updated for each new frame in the middle loop, without a new selection being created or destroyed. There is no equivalent option for an angular "within" so angular histogramming occurs more traditionally via a loop over lipids.  

proc loop_over_shells {species headname tailname lipidbeads_selstr low_f upp_f low_f_avg upp_f_avg} {
    global params
    set delta_frame [expr ($params(end_frame) - $params(start_frame)) / $params(dt)]
    for {set ri $params(Rmin)} { $ri<=$params(Rmax)} { set ri [expr $ri + $params(dr)]} {
        #loop over shells
        puts "Now on shell {$ri [expr ${ri}+$params(dr)]}"
        set rf [expr $ri + $params(dr)]
        set rf2 [expr $rf*$rf]
        set ri2 [expr $ri*$ri]
        set shell [atomselect top "(resname $species) and ((x*x + y*y < $rf2) and  (x*x + y*y > $ri2)) and $lipidbeads_selstr"]
        #puts [$shell num]		
        set theta_bin [loop_over_frames $shell "resname $species" $headname $tailname $lipidbeads_selstr $params(start_frame) $params(end_frame)  $ri $rf $low_f $upp_f]
        set theta_bin_high [lindex $theta_bin 1]
        set theta_bin_low [lindex $theta_bin 0]
        $shell delete	
        set time_avg_upper [vecscale $theta_bin_high [expr 1.0 / (1.0 * $delta_frame)]]
        set time_avg_lower [vecscale $theta_bin_low [expr 1.0 / (1.0 * $delta_frame)]]
        binningTools::output_bins $upp_f_avg $ri $rf "$time_avg_upper" 
        binningTools::output_bins $low_f_avg $ri $rf "$time_avg_lower" 
    }
}


proc set_parameters { config_file_script } {
    global params
    array unset params
    global params
    array set params {
        leaflet_sorting_algorithm 2
        center_and_align 0
        use_qwrap 0
        utils ./helpers 
        dt 1
        leaflet_reassign_interval 1
        start_frame 0  
        backbone_selstr "name BB" 
        protein_selstr "name BB SC1 to SC4"
        lipids {"POPG"}
        headnames {"PO4"}
        tailnames {"C4"}
        lipidbeads_selstrs {"all"}
        acylchain_selstrs {"all"}        
        chainlist {A B C D E}
        helixlist {1 2 3 4}
        helix_assignment_script "./helpers/assign_helices_2BG9_CG_lms2.tcl" 
        midplane_selstr "occupancy 1 to 4"        
        Rmax 20.
        Rmin 0.
        dr 5 
        Ntheta 50        
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


;#The main function that initializes, constructs the densities for each lipid species, and outputs to file. 

proc polarDensityBin { config_file_script } { 
    ;#read parameters
    #source $config_file_script
    global params
    set_parameters $config_file_script
    source $params(utils)/binning_tools.tcl
    source $params(utils)/math_tools.tcl
    source $params(utils)/orientation_membOrganization.tcl
    if {$params(use_qwrap) == 1} {load $params(utils)/qwrap.so}
    set backbone_selstr $params(backbone_selstr) ;#only necessary for backwards compatibility 
    set protein_selstr $params(protein_selstr) ;#only necessary for backwards compatibility 
    source $params(helix_assignment_script)
    foreach species $params(lipids) lipidbeads_selstr $params(lipidbeads_selstrs) acylchain_selstr $params(acylchain_selstrs) headname $params(headnames) tailname $params(tailnames) {
        set outfile "$species"
        set sel [atomselect top "resname $species"]
        set sel_num [$sel num]
        
        if {$sel_num == 0} {
            error "No lipid of species $species"
        }
        ;# Center's system (weak hack)
        if {$params(center_and_align) == 1} {
            membOrganization::center_and_wrap_system "occupancy $params(helixlist) and $params(backbone_selstr)"
            membOrganization::center_and_wrap_system "occupancy $params(helixlist) and $params(backbone_selstr)"
            membOrganization::center_and_wrap_system "occupancy $params(helixlist) and $params(backbone_selstr)"
            ;# aligns protein
            membOrganization::Align "occupancy $params(helixlist) and $params(backbone_selstr)"
        }
        ;# outputs protein positions
        output_helix_centers
        ;# initialize some constants
        set area [mathTools::get_avg_area top]
        set nframes [molinfo top get numframes]
        if { $params(start_frame) > $nframes } {
            puts "Error: specified start frame $params(start_frame) is greater than number of frames $nframes" 
            set start_frame $nframes
        }
        if { $params(end_frame) > $nframes } {
            puts "Warning: specified end frame $params(end_frame) is greater than number of frames; setting end frame to $nframes" 
            set end_frame $nframes
        }
        $sel delete

        puts "Acyl Chain:\t$species"
        set low_f [open "${outfile}.low.dat" w]
        set upp_f [open "${outfile}.upp.dat" w]
        set low_f_avg [open "${outfile}.low.avg.dat" w]
        set upp_f_avg [open "${outfile}.upp.avg.dat" w]

        set totals [lipid_assessment::frame_leaflet_assignment "resname $species" $headname $tailname $lipidbeads_selstr $params(start_frame) $params(start_frame)]
        
        foreach lu [list $low_f $upp_f] avgfile [list $low_f_avg $upp_f_avg] leaf_total $totals {
            set leaflet_str [lindex $leaf_total 0]
            set expected_beads [lindex $leaf_total 1]
            set expected_lipids [lindex $leaf_total 2]
            set expected_bead_density [expr 1.0 * $expected_beads/$area]
            puts "#Lipid species $species in $leaflet_str leaflet: ${expected_lipids} molecules, Num beads : ${expected_beads} beads,  Average Area : [format {%0.0f} $area] A^2, Expected Bead Density : [format {%0.5f} [expr $expected_bead_density]]/A^2, Average Chain : [membOrganization::avg_acyl_chain_len "resname $species" $acylchain_selstr] beads, dr*dtheta : [format {%0.5f} [expr $params(dr)*[mathTools::DtoR $params(dtheta)]]] "
            puts $lu "#Lipid species $species in $leaflet_str leaflet: ${expected_lipids} molecules, Num beads : ${expected_beads} beads,  Average Area : [format {%0.0f} $area] A^2, Expected Bead Density : [format {%0.5f} [expr $expected_bead_density]]/A^2, Average Chain : [membOrganization::avg_acyl_chain_len "resname $species" $acylchain_selstr] beads, dr*dtheta : [format {%0.5f} [expr $params(dr)*[mathTools::DtoR $params(dtheta)]]] "
            puts $avgfile "#Lipid species $species in $leaflet_str leaflet: ${expected_lipids} molecules, Num beads : ${expected_beads} beads,  Average Area : [format {%0.0f} $area] A^2, Expected Bead Density : [format {%0.5f} [expr $expected_bead_density]]/A^2, Average Chain : [membOrganization::avg_acyl_chain_len "resname $species" $acylchain_selstr] beads, dr*dtheta : [format {%0.5f} [expr $params(dr)*[mathTools::DtoR $params(dtheta)]]] "
        }
        puts "Processing frames, starting at frame $params(start_frame) and ending at frame $params(end_frame)."
        lipid_assessment::trajectory_leaflet_assignment "resname $species" $headname $tailname $lipidbeads_selstr         
        ;#the core calculation 
        loop_over_shells $species $headname $tailname $lipidbeads_selstr $low_f $upp_f $low_f_avg $upp_f_avg  
        
        close $low_f
        close $upp_f
        close $low_f_avg
        close $upp_f_avg
    }
}    