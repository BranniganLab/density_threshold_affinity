set selall [atomselect top {all}]
$selall set occupancy 0.0
$selall set user 0.0
$selall set chain "X"
$selall delete

set helix_code_list [list 1 2 3 4]
set proBeads [atomselect top "name BB SC1 to SC4"]
set nProtein [expr [$proBeads num]/3892]
$proBeads delete


set M1_start    201
set M1_end	    226
set M2_start    228
set M2_end	250
set M3_start	260
set M3_end	285
set M4_start	293
set M4_end	320

;# Start and end residue numbers for M1 to M4
set M1_start_list [[atomselect top "resid $M1_start and name BB"] get residue]
set M4_end_list   [[atomselect top "resid $M4_end and name BB"] get residue]

puts $M1_start_list
puts $M4_end_list
foreach start $M1_start_list end $M4_end_list chainId [list E D C B A] {
    set tempsel [atomselect top "residue $start to $end"]
    $tempsel set chain $chainId
    $tempsel delete
}


[atomselect top "resid $M1_start to $M1_end and name BB SC1 to SC4"] set occupancy 1
[atomselect top "resid $M2_start to $M2_end and name BB SC1 to SC4"] set occupancy 2
[atomselect top "resid $M3_start to $M3_end and name BB SC1 to SC4"] set occupancy 3
[atomselect top "resid $M4_start to $M4_end and name BB SC1 to SC4"] set occupancy 4
