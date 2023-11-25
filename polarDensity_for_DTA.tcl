
 package require pbctools
 set UTILS "/u2/home_u2/lms464/github/JPC_Special/common/utils" 
 set QWRAP "/u1/home/lms464/lms464/github/qwrap"
# 
source $UTILS/BinTools.tcl
# ;# TODO What is the point of outputing each lipid species to a different file?
# 
load ${QWRAP}/qwrap.so
# 
# ;#Lipid_Saturation_HeadG are a series of macros to parse Martini lipids
# source ${UTILS}/Lipid_Saturation_HeadG.tcl


#Grace Brannigan 7/2018

#Sample Use:
#Assuming trajectory of interest is loaded AND top.  All frames will be used so any frames to be ignored would need to be unloaded
#set HeadNames [atomselect top "name PO4 ROH B2"] ;#B2 is DDM
#set lipids [lsort -unique [$HeadNames get resname]]
#$HeadNames delete
#set RMax 40.
#set RMin 5.
#set dr 2.
#set Ntheta 30
#foreach lip in $lipids {
#polarDensityBin $lip.dat $lip $Rmin $Rmax $dr $Ntheta	
#}

###############################################Sample plotting within python
#
#Ntheta = 30
#data = np.loadtxt('DPPC.dat',skiprows=2)
#rad = data[:,1] + (data[:,1]-data[:,0])/2.0
#the = np.linspace(0,2*np.pi,Ntheta +1)
#theta,radius=np.meshgrid(the,rad)
#density = data[:,3:]/radius 
#plt.figure(figsize = (5,5))
#plt.subplot(projection="polar")
#plt.pcolormesh(theta,radius,density,cmap="RdBu",zorder=0,edgecolors='k',lw=.001)
#plt.show()

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

# Sum_list
#
# Sums over a list
# Arguments:
#   list: a list of numbers
# Outputs:
#   float: the sum over the entire list
#
# Issues:
#    Function name violates style guide
proc Sum_list {list_in} {
    set list_out 0
    foreach li $list_in {
        set list_out [expr 1.0*$list_out+$li]
    }
    return $list_out
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
# Issues:
#    Headgroup beads should be an argument
proc z_mid {init_frm nframes} {
    set z_list {}
    for {set frm ${init_frm}} {${frm} < ${nframes}} {incr frm} {
        set mid [atomselect top "name PO4 ROH C3 PO41" frame $frm]
        lappend z_list [lindex [measure center $mid weight mass] 2]
        $mid delete
    }
    return [expr 1.0*[vecsum $z_list]/([llength $z_list]) ]
}

### Debug helpers ###
# Shell_Test
#
# Tests that the sum over theta is the same as the sum over r. (I think - ES)
# Arguments:
#   float: the sum over r (?)
#   list of floats: counts for theta bins
# Result:
#   Prints a string indicating whether or not the counts are consistent
#
# Issues:
#    Should use actual error handling
proc Shell_Test {shell_count theta_bin_counts } {
    set theta_bin_total [Sum_list [lindex $theta_bin_counts 1]] 
    if {$shell_count == $theta_bin_total} {
        puts "Counting appears to be consistent"
    } else {
        puts "Counting appears to be inconsistent.."
    }

}

# Sum_Shell_Test
#
# Tests that the sum over theta is the same as the sum over r. (I think - ES)
# Arguments:
#   str: lipid species to analyze - seltext
#   float: Rmin, inner radius of shell
#   float: Rmax, outer radius of shell
#   float: dr, distance between Rmin and Rmax
#   int: frm, the frame to analyze
#   int: sel_num, the expected number of beads in the shell
# Result:
#   Prints the difference between the selection and the sum
#
# Issues:
#    Should use actual error handling
#    Should be made into a proper test
proc Sum_Shell_Warning {species Rmin Rmax dr frm sel_num} {

	set total_beads 0 
	for {set ri $Rmin} { $ri<=${Rmax}} {set ri [expr $ri + $dr]} {
		set rf [expr $ri + $dr]
        set rf2 [expr $rf*$rf]
        set ri2 [expr $ri*$ri]
		set test_shell [atomselect top "($species) and ((x*x + y*y < $rf2) and  (x*x + y*y > $ri2)) " frame $frm]
		incr total_beads [$test_shell num]
	} 
	
	puts "[expr abs($sel_num-$total_beads)]"
}

# Species_Total_Warning
#
# Makes sure not all lipids are in one shell
# Arguments:
#   int: The number of lipids/beads in the selection
#   int: The number of lipids/beads in the shell
# Result:
#   Prints a warning message if the two ints are the same.
#
# Issues:
#    Should use actual error handling
#    Should be made into a proper test
proc Species_Total_Warning {sel_num shell_count} {
    if {$shell_count == $sel_num} {
        puts "Warning: One shell appears to contain all lipids!"
        puts "Not exceptable! Exiting."
        exit
    }
}

# Real_vs_Expected
#
# Compares a bin's density with the expected density. Prints a diagnostic string.
# Arguments:
#   float: the expected density
#   float: bin_counts, the counts from the bins of the shell
#   float: ri, the inner radius of the shell
#   float: rf, the outer radius of the shell
# Result:
#   Prints a diagnostic string indicating if the shell is enriched, depeleted, or randomly mixed
#
# Issues:
#    It's not clear if this is a sanity check (sniff test) or useful data. -ES
proc Real_vs_Expected {expected bin_counts ri rf} {
    set scalar [expr 1.0/(3.14159265358979323846*($rf**2-$ri**2))]
    set test_list [vecscale $bin_counts $scalar]
    set avg_test_list [expr 1.0*[Sum_list $test_list]/[llength $bin_counts]]
    set ratio [expr 1.0*$avg_test_list/$expected]
    if {$ratio < .8} {
        puts "This shell appears to be depleted\nAvg: $avg_test_list\tExpected: $expected\tRatio: $ratio"
    } elseif {$ratio > 1.1} {
        puts "This shell apears to be enriched\nAvg: $avg_test_list\tExpected: $expected\tRatio: $ratio"
    } else {
        puts "Shell appears to be randomly mixed"
    }
}

# Determines a specific lipids leaflet



# Ouputs position of the centered protein in a membrane
# accross both leaflets
proc Protein_Position {{a ""}} {
    set chain_names [list "A" "B" "C" "D" "E"]
    set zed [z_mid 0 20]
	set occupancy [list 1 2 3 4]
	foreach eq {"<" ">"} eqtxt {"lwr" "upr"} {
		set fout [open "/u2/home_u2/lms464/github/JPC_Special/tasks/17_Aim1/Data/Protein${a}_coords_${eqtxt}.dat" w]
        puts $fout  "# chain A ooc1r occ1the occ2r occ2the... "
        foreach chnm $chain_names {
            foreach occ $occupancy {
                set sel [atomselect top "(chain ${chnm}) and (occupancy $occ and name BB) and (z ${eq} $zed)" frame 0]
                set com [measure center $sel weight mass]
                $sel delete
                set x [lindex $com 0]
                set y [lindex $com 1]
                set r [expr sqrt($x*$x+$y*$y)]
                set theta [get_theta $x $y]
                puts "chain ${chnm} and occupancy $occ $r $theta"

                puts -nonewline $fout "$r $theta "
            }
            puts $fout ""
        }
        close $fout
    }
}

proc avg_acyl_chain_len {species} {
    
    set acyl_num 0
    set sel [atomselect top "$species"]
    set sel_resname [lsort -unique [$sel get resname]]
    #set sel_num [llength [lsort -unique [$sel get resname]]]
    $sel delete
    foreach res $sel_resname {
        set sel [atomselect top "${species} and (resname $res) and (not name NH3 NC3 GL1 GL2 AM1 AM2 PO4 CNO CN0 C1 C2 C3)"]
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

proc Center_System {inpt} {
    puts "${inpt}"
    # confirms your box is either square or paraelleogram-ish
    # will preform qwrap or pbc wrap depending

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
        
        if {([lindex $pbc_angles 0]!=90.0) && ([lindex $pbc_angles 0]!=90.0) && ([lindex $pbc_angles 0]!=90.0)} {
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

proc resnamer {input} {

    
    # adds resname if the input is DPPC, CHOL, PUPI...
    
    set out ""
    if {[string length $input] == 4 && $input != "chol"} { 
        set out "resname $input"
    } else {
        set out "$input"
    }
    return $out
}

proc output_bins {fl  ri rf dtheta bins} {
    puts -nonewline $fl "[format {%0.2f} $ri]  [format {%0.2f} $rf] [format {%0.2f} $dtheta]  " 
    puts $fl "$bins" 
}

# LMS keep here or a new file (please move it to a new file!)
# proc bin_over_frames {shell species dtheta sample_frame nframes dt } {
#     set theta_bin_high [list ]
#     set theta_bin_low [list]
#     set shel_count 0
#     for {set frm $sample_frame} {$frm < ${nframes}} {incr frm $dt} {
#     #loop over frames
#         $shell frame $frm
#         $shell update 
#         set indexs [$shell get index]
        
#         ;# trying to implement meaningful leaflet seperation
#         ;# change 1
#         set resids [$shell get resid]
#         ;# end of new change 1
        
#         set nShell [$shell num]
#         ;#Species_Total_Warning $sel_num $nShell
#         set shel_count [expr $shel_count + $nShell]
#         set ti_list {}
#         ;#new_var change 3

#         ;# change 2 adding to the foreach command
#         foreach indx $indexs resd $resids {
#             set high_low 0
#             #loop over lipids in the shell
#             set a "$species and index $indx"
#             set b "(resid $resd) and (name PO4 ROH)"
#             ;# change 4
#             set thislipid [atomselect top $a frame $frm]
#                 if {[$thislipid get name] == "PO4"} {
#                  continue
#                 }
#             ;# change 5
#             set high_low [local_mid_plane $b $resd "PO4 ROH" $frm]
            
#             set x [$thislipid get x]
#             set y [$thislipid get y]
#             $thislipid delete
#             set theta [get_theta $x $y]
#             set ti [expr int($theta/$dtheta)] 
#             #determine theta bin
#             if {$high_low > 0} {
#                 lappend theta_bin_low $ti
#             } else {
#                 lappend theta_bin_high $ti
#             }
#             ;#lappend theta_bin_list $ti
#             #add an instance of this theta bin to a list for later counting
#         }           
#     }
#     return [list ${theta_bin_high} ${theta_bin_low} $shel_count]

# }

proc bin_over_frames {shell species dtheta sample_frame nframes Ntheta dt ri rf fupper flower} {
    set theta_bin_high [lrepeat [expr $Ntheta+1] 0]
    set theta_bin_low [lrepeat [expr $Ntheta+1] 0]
    for {set frm $sample_frame} {$frm < ${nframes}} {incr frm $dt} {
        #loop over frames
	    #puts $frm
        $shell frame $frm
        $shell update 
        set singleFrame_counts [bin_frame $shell $species $dtheta $frm ]
        # you'll need to create bin_frame, using lines 284-325 (or around those) of your previous code
        set singleFrame_upper [lindex $singleFrame_counts 0] 
	    #puts $singleFrame_upper
        #I assume here that bin_frame returns upper and lower as two lists inside another list, you can do it however
        set singleFrame_lower [lindex $singleFrame_counts 1]
        set theta_bins [theta_histogram $singleFrame_upper $singleFrame_lower $Ntheta]

        # should be fixed, do not change [lrepeat [expr $Ntheta+1] to [lrepeat [expr $Ntheta] 
        if { [llength $theta_bin_high] != [llength [lindex $theta_bins 0]] } {
            error "theta_bin_high/low and theta_bins do not have the same length."
        }
        set theta_bin_high [vecadd $theta_bin_high [lindex $theta_bins 0] ]
        #puts [lindex $theta_bins 1]
	    set theta_bin_low [vecadd $theta_bin_low [lindex $theta_bins 1]]
	    #puts $theta_bin_low
        #TODO MAKE A SWITCH
	#output_bins $fupper $ri $rf $dtheta [lindex $theta_bins 0] 
        #open fupper before the loop starts and close afterwards
        #output_bins $flower $ri $rf $dtheta [lindex $theta_bins 1] 
        #same thing     
    }
  return [list ${theta_bin_high} ${theta_bin_low}]
}

proc local_mid_plane {atsel_in frame_i} {
    set temp_sel [atomselect top "(name PO4 ROH) and (pbwithin 50 of $atsel_in) and not ($atsel_in)" frame $frame_i]
    set mid_point [lindex [measure center $temp_sel weight mass] 2]
    $temp_sel delete
    set sel_resid [atomselect top "$atsel_in" frame $frame_i]
    set resid_z [lindex [lindex [${sel_resid} get {x y z}] 0] 2]
    $sel_resid delete
        
    if {$mid_point < $resid_z} {
        return 0
    } else {
        return 1
    }
}

# proc local_mid_plane {atsel_in frame_i} {

# 	;# Method ONLY works under the assumption lipids are hetero-acidic
# 	;# having a saturate or mono-unsaturated  

#     set sel_resid [atomselect top "$atsel_in" frame $frame_i]
#     set ind 1
#     ;# Says choose the PO4 bead. PA dosen't have a 2 bead head group, only PO4
#     if { [string range [lsort -unique [$sel_resid get resname]] end-1 end] == "PA" } {
#     	set ind 0
#     }
#     set sel_Z [${sel_resid} get z] 
#     ;#set sel_com [measure center ${sel_resid} weight mass]
#     $sel_resid delete
# 	if {[lindex ${sel_Z} $ind] < [lindex ${sel_Z} end] } { 
# 		return 1 
# 	} else { 
# 		return 0
# 	}
# }



# does what it says it does, bins over a single frame
proc bin_frame {shell species dtheta frm } {
    set indexs [$shell get index]
    set resids [$shell get resid]
    set nShell [$shell num]
    set theta_high_out [list]
    set theta_low_out [list]
    set resd_old 0
    set high_low 0
    #set shel_count [expr $shel_count + $nShell]
    foreach indx $indexs resd $resids {
        #loop over lipids in the shell
        set a "($species and index $indx)"
        set b "(resid $resd)" ;#and (name PO4 ROH)
        # change 4
        set thislipid [atomselect top $a frame $frm]
        if {[string length ${species}] == 2} {
	    	if {[$thislipid get name] == "PO4"} {
	        	continue
	    	}
    	}
        # change 5

        if {${resd_old} != ${resd}} {
        	set high_low [local_mid_plane $b  $frm]
        }
        set x [$thislipid get x]
        set y [$thislipid get y]
        $thislipid delete
        set theta [get_theta $x $y]
        set ti [expr int($theta/$dtheta)] 
        #determine theta bin
        if {$high_low > 0} {
            lappend theta_low_out $ti
        } else {
            lappend theta_high_out $ti
        }

    }
    
    return [list $theta_high_out $theta_low_out] 
}

# FAR more useful than the other version (theta clean up)
proc theta_histogram {singleFrame_upper singleFrame_lower Ntheta } {
    
    set theta_bin_out [list]

    foreach ud [list $singleFrame_upper $singleFrame_lower] {
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

# TODO I don't think I need this function anymore
proc theta_clean_up {theta_bin_high theta_bin_low shel_count  Ntheta delta_frame low_f upp_f} {
    
    theta_bin_out [list ]

    foreach ud [list $theta_bin_low $theta_bin_high] {
        #Species_Total_Warning $sel_num $shel_count
        puts "Cleaning up for shell $ri to $rf"
        #cleanup and output 
        set theta_bin_counts [lcount $ud]
        #Shell_Test $shel_count $theta_bin_counts
        set theta_bin_time_averages {}
        for {set ti 0} { $ti<=$Ntheta} {incr ti 1} {
            set tindex [lsearch [lindex $theta_bin_counts 0] $ti]
            if { $tindex >= 0} {
                set time_average [expr 1.0 * [lindex [lindex $theta_bin_counts 1] $tindex]/(1.0*($delta_frame))] 
            } else { 
                set time_average 0.0
            }
            lappend theta_bin_time_averages $time_average
        }
        lappend theta_bin_out $theta_bin_time_averages
    }
    return $theta_bin_time_averages
}

### polarDensity Funciton ###

proc polarDensityBin { outfile species Rmin Rmax dr Ntheta} {
    #if {$species == "CHOL"} { set species "resname CHOL" }
	#source /u2/home_u2/lms464/github/JPC_Special/common/grace/assign_helices_3RQW_CG_lms.tcl;#assign_helices_2BG9_CG_lms2.tcl
	source /u2/home_u2/lms464/github/JPC_Special/common/grace/assign_helices_2BG9_CG_lms2.tcl
    
    #set Rmin 0
    #set Rmax 36

    set species [resnamer ${species}]
    
    set sel [atomselect top "$species"]
	set sel_num [$sel num]
	
	if {$sel_num == 0} {
        error "No lipid saturation set exists"
	}
	
	;# TODO CHANGE THIS BACK
 	Center_System "occupancy 1 to 4 and name BB"
    Center_System "occupancy 1 to 4 and name BB"
    Center_System "occupancy 1 to 4 and name BB"

 	Align "occupancy 1 to 4 and name BB"
    Protein_Position   
	set dt 1
    set area [get_avg_area top]
	set nframes [molinfo top get numframes]
    $sel delete
	puts "Acyl Chain:\t$species"
	set low_f [open "${outfile}.low.dat" w]
    set upp_f [open "${outfile}.upp.dat" w]
	set dtheta [expr 360.0/(1.0*($Ntheta))]
    #Center_System "name PO4"
	
    foreach lu [list $low_f $upp_f] zed [list "(z<0)" "(z>0)"] {
        set sel [ atomselect top "(($species) and $zed) and (name PO4 ROH)"  frame 0]
        set sel_num [llength [lsort -unique [$sel get resid] ] ]
        if {$sel_num < 1} {
            set num_beads 0
        	set expected 0
        } else {
	        set sel_resid [lsort -unique [$sel get resid] ]
	        $sel delete
	        set beads [atomselect top "${species} and (resid $sel_resid) and (not name PO4)" frame 0]
	        set num_beads [$beads num]
	        set expected [expr 1.0 * $num_beads/$area]
	        $beads delete
    	}
        puts "#Lipid species $species : ${sel_num} molecules, Num beads : ${num_beads} beads,  Average Area : [format {%0.0f} $area] A^2, Expected Density : [format {%0.5f} [expr $expected]]/A^2, Average Chain : [avg_acyl_chain_len ${species}] beads, dr*dtheta : [format {%0.5f} [expr $dr*[DtoR $dtheta]]] "
	    puts $lu "#Lipid species $species : ${sel_num} molecules, Num beads : ${num_beads} beads,  Average Area : [format {%0.0f} $area] A^2, Expected Density : [format {%0.5f} [expr $expected]]/A^2, Average Chain : [avg_acyl_chain_len ${species}] beads, dr*dtheta : [format {%0.5f} [expr $dr*[DtoR $dtheta]]] "
    }
    #Center_System "occupancy 1 to 4 and name BB"
    #Align "occupancy 1 to 4 and name BB"

    
	#unset lipsize
	set sample_frame 100
	set delta_frame [expr ($nframes - $sample_frame) / $dt]
	for {set ri $Rmin} { $ri<=${Rmax}} { set ri [expr $ri + $dr]} {
		#loop over shells
		puts "Ring {$ri [expr ${ri}+${dr}]}"
		set rf [expr $ri + $dr]
		set rf2 [expr $rf*$rf]
		set ri2 [expr $ri*$ri]
		set shell [atomselect top "($species) and ((x*x + y*y < $rf2) and  (x*x + y*y > $ri2)) and (not name PO4)"]
		#selects lipids in the radial shell
		#set shel_count 0
        #set theta_bin_high {}
        #set theta_bin_low {}		
        set theta_bin [bin_over_frames $shell $species $dtheta $sample_frame $nframes $Ntheta $dt $ri $rf $upp_f $low_f]
        set theta_bin_high [lindex $theta_bin 0]
        set theta_bin_low [lindex $theta_bin 1]
        #puts ${theta_bin_high}
        #set shel_count [expr $shel_count + [lindex $theta_bin 2]]
        $shell delete	
        puts ""
        #theta_bin_averages [theta_clean_up $theta_bin_high $theta_bin_low $shel_count $Ntheta $dtheta $delta_frame]
        set time_avg_upper [vecscale $theta_bin_high [expr 1.0 / (1.0 * $delta_frame)]]
        set time_avg_lower [vecscale $theta_bin_low [expr 1.0 / (1.0 * $delta_frame)]]
		puts ""
        #foreach tbu [lindex $theta_bin_averages 0] tbl [lindex $theta_bin_averages 1] ;#{
        output_bins $upp_f $ri $rf $dtheta "$time_avg_upper" 
        output_bins $low_f $ri $rf $dtheta "$time_avg_lower" 
        #}
	}
	close $low_f
	close $upp_f
}
