set selall [atomselect top {all}]
$selall set occupancy 0.0
$selall set user 0.0
$selall delete

set helix_code_list [list 1 2 3 4 5 6 7]
set proBeads [atomselect top "name BB SC1 to SC4"]
set nProtein [expr [$proBeads num]/3892]
$proBeads delete


set TM1_start   1 
set TM1_end	    31
set TM2_start   36
set TM2_end	    65
set TM3_start	70
set TM3_end	    105
set TM4_start	115
set TM4_end	    140
set TM5_start	165
set TM5_end	    200
set TM6_start	205
set TM6_end	    235
set TM7_start	245
set TM7_end	    265

[atomselect top "resid $TM1_start to $TM1_end and name BB SC1 to SC4"] set occupancy 1
[atomselect top "resid $TM2_start to $TM2_end and name BB SC1 to SC4"] set occupancy 2
[atomselect top "resid $TM3_start to $TM3_end and name BB SC1 to SC4"] set occupancy 3
[atomselect top "resid $TM4_start to $TM4_end and name BB SC1 to SC4"] set occupancy 4
[atomselect top "resid $TM5_start to $TM5_end and name BB SC1 to SC4"] set occupancy 5
[atomselect top "resid $TM6_start to $TM6_end and name BB SC1 to SC4"] set occupancy 6
[atomselect top "resid $TM7_start to $TM7_end and name BB SC1 to SC4"] set occupancy 7

[atomselect top "resid $TM1_start to $TM7_end and name BB SC1 to SC4"] set chain A 