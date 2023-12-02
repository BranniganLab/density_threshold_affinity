set selall [atomselect top {all}]
$selall set occupancy 0.0
$selall set user 0.0
$selall delete

set helix_code_list [list 1 2 3 4]
set proBeads [atomselect top "name BB SC1 to SC4"]
set nProtein [expr [$proBeads num]/3892]
$proBeads delete

;# Start and end resid numbers for M1 to M4
set M1_start_list [list 211 208 210 211 212]
set M1_end_list   [list 236 233 236 236 236]
set M2_start_list [list 243 241 242 243 244]
set M2_end_list   [list 270 267 269 270 271]
set M3_start_list [list 274 274 274 274 274]
set M3_end_list   [list 302 302 302 302 302]
set M4_start_list [list 308 307 307 308 307]
set M4_end_list   [list 700 336 700 700 700]

;# lists of each occupancy starting ending points
set res_start_list_2D [list $M1_start_list $M2_start_list $M3_start_list $M4_start_list] 
set res_end_list_2D [list $M1_end_list $M2_end_list $M3_end_list $M4_end_list] 

set A_num 787
set B_num 774
set C_num 773
set D_num 787
set E_num 771
set hold1 1
set hold2 0

set serial_start_list [list ]
set serial_end_list [list ]

for {set i 0} {$i < $nProtein} {incr i} {
    foreach num [list $A_num $B_num $C_num $D_num $E_num] {
        set hold2 [expr $hold2 + $num]
        lappend serial_start_list $hold1
        lappend serial_end_list $hold2
        set hold1 [expr $num + $hold1]
    }
}   

set chainlist [list A B C D E F G H I J K L M N O P Q R S T U V W "1" Y Z a b c d e f g h i j k l m n o p q r s t u v w "2" "3" "4" "5" "6" "7" "8" "9" "0"]

;# Names chains
for {set i 0} {$i < [expr 5*$nProtein]} {incr i} {
  set serialstart [lindex $serial_start_list $i]
  set serialend [lindex $serial_end_list $i]
  set chainsel [atomselect top "serial $serialstart to $serialend"]
  $chainsel set chain [lindex $chainlist $i]
  $chainsel delete
}
;# sets up occupancy
for {set lranger 0} {$lranger <= [expr 5*$nProtein - 1]} {incr lranger 5} {
    set serial_start_block_5 [lrange $serial_start_list $lranger [expr 4 + $lranger]]
    set serial_end_block_5 [lrange $serial_end_list $lranger [expr 4 + $lranger]]
    foreach helix_code $helix_code_list res_start_list $res_start_list_2D res_end_list $res_end_list_2D {
    ;#    puts "\nDoing helix $helix_code"
        set subunit 0
        foreach serial_start $serial_start_block_5 serial_end $serial_end_block_5 res_start $res_start_list res_end $res_end_list {
        
           ;# puts "$serial_start $serial_end | $res_start $res_end"
            set tmp_sel [atomselect top "name BB SC1 SC2 SC3 SC4 and resid $res_start to $res_end"]
            set tmp_serials [$tmp_sel get serial]
            $tmp_sel delete
            set serials ""   
            foreach serial $tmp_serials   {
                if {($serial < $serial_end) && ($serial > $serial_start) } {
                    lappend serials $serial
                }
            } 
            set sel [atomselect top "same residue as (serial $serials and resid $res_start to $res_end)"]
            $sel set user $helix_code
            $sel set occupancy $helix_code
            incr subunit
            $sel delete
        }
    }
    unset res_start res_end
}
