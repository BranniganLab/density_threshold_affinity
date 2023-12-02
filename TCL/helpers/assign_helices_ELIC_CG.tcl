set selall [atomselect top {all}]
$selall set occupancy 0.0
$selall set user 0.0
$selall delete

set helix_code_list [list 1 2 3 4]
set proBeads [atomselect top "name BB SC1 to SC4"]
set nProtein [expr [$proBeads num]/3892]
$proBeads delete


set M1_start    192 
set M1_end	    216
set M2_start    218
set M2_end	243
set M3_start	251
set M3_end	272
set M4_start	286
set M4_end	307

;# Start and end residue numbers for M1 to M4
set M1_start_list [[atomselect top "resid $M1_start and name BB"] get residue]
set M1_end_list   [[atomselect top "resid $M1_end and name BB"] get residue]
set M2_start_list [[atomselect top "resid $M2_start and name BB"] get residue]
set M2_end_list   [[atomselect top "resid $M2_end and name BB"] get residue]
set M3_start_list [[atomselect top "resid $M3_start and name BB"] get residue]
set M3_end_list   [[atomselect top "resid $M3_end and name BB"] get residue]
set M4_start_list [[atomselect top "resid $M4_start and name BB"] get residue]
set M4_end_list   [[atomselect top "resid $M4_end and name BB"] get residue]

puts $M1_start_list
puts $M4_end_list
foreach start $M1_start_list end $M4_end_list chainId [list A B C D E] {
    set tempsel [atomselect top "residue $start to $end"]
    $tempsel set chain $chainId
}


[atomselect top "resid $M1_start to $M1_end and name BB SC1 to SC4"] set occupancy 1
[atomselect top "resid $M2_start to $M2_end and name BB SC1 to SC4"] set occupancy 2
[atomselect top "resid $M3_start to $M3_end and name BB SC1 to SC4"] set occupancy 3
[atomselect top "resid $M4_start to $M4_end and name BB SC1 to SC4"] set occupancy 4




;# lists of each occupancy starting ending points
# set res_start_list_2D [list $M1_start_list $M2_start_list $M3_start_list $M4_start_list] 
# set res_end_list_2D [list $M1_end_list $M2_end_list $M3_end_list $M4_end_list] 

# set A_num 787
# set B_num 774
# set C_num 773
# set D_num 787
# set E_num 771
# set hold1 1
# set hold2 0

# set serial_start_list [list ]
# set serial_end_list [list ]

# for {set i 0} {$i < $nProtein} {incr i} {
#     foreach num [list $A_num $B_num $C_num $D_num $E_num] {
#         set hold2 [expr $hold2 + $num]
#         lappend serial_start_list $hold1
#         lappend serial_end_list $hold2
#         set hold1 [expr $num + $hold1]
#     }
# }   

# set chainlist [list A B C D E F G H I J K L M N O P Q R S T U V W "1" Y Z a b c d e f g h i j k l m n o p q r s t u v w "2" "3" "4" "5" "6" "7" "8" "9" "0"]

# ;# Names chains
# for {set i 0} {$i < [expr 5*$nProtein]} {incr i} {
#   set serialstart [lindex $serial_start_list $i]
#   set serialend [lindex $serial_end_list $i]
#   set chainsel [atomselect top "serial $serialstart to $serialend"]
#   $chainsel set chain [lindex $chainlist $i]
#   $chainsel delete
# }
# ;# sets up occupancy
# for {set lranger 0} {$lranger <= [expr 5*$nProtein - 1]} {incr lranger 5} {
#     set serial_start_block_5 [lrange $serial_start_list $lranger [expr 4 + $lranger]]
#     set serial_end_block_5 [lrange $serial_end_list $lranger [expr 4 + $lranger]]
#     foreach helix_code $helix_code_list res_start_list $res_start_list_2D res_end_list $res_end_list_2D {
#     ;#    puts "\nDoing helix $helix_code"
#         set subunit 0
#         foreach serial_start $serial_start_block_5 serial_end $serial_end_block_5 res_start $res_start_list res_end $res_end_list {
        
#            ;# puts "$serial_start $serial_end | $res_start $res_end"
#             set tmp_sel [atomselect top "name BB SC1 SC2 SC3 SC4 and resid $res_start to $res_end"]
#             set tmp_serials [$tmp_sel get serial]
#             $tmp_sel delete
#             set serials ""   
#             foreach serial $tmp_serials   {
#                 if {($serial < $serial_end) && ($serial > $serial_start) } {
#                     lappend serials $serial
#                 }
#             } 
#             set sel [atomselect top "same residue as (serial $serials and resid $res_start to $res_end)"]
#             $sel set user $helix_code
#             $sel set occupancy $helix_code
#             incr subunit
#             $sel delete
#         }
#     }
#     unset res_start res_end
# }
