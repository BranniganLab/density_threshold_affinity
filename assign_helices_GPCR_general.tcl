set helix_code_list [list 1 2 3 4 5 6 7]

#canonical numbering 
set start_list [list 32 67 101 146 196 267 305]
set end_list [list 62 96 136 171 231 300 329]

#gromacs numbering, comment out for canonical numbering
set start_list [list 1 36 70 115 165 205 245]
set end_list [list 31 65 105 140 200 235 265]


set res_sel_str "alpha" ;#atomistic
set res_sel_str "name BB SC1 to SC4" ;#Martini, comment out for atomistic 


set sel [atomselect top "all"]
$sel set occupancy 0

foreach start_res $start_list end_res $end_list helix_code $helix_code_list {
	puts "$start_res $end_res $helix_code"
	[atomselect top "resid $start_res to $end_res and $res_sel_str"] set occupancy $helix_code
}

[atomselect top "$res_sel_str"] set chain A 