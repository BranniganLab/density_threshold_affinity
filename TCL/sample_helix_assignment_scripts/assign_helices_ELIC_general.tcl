set selall [atomselect top {all}]
$selall set occupancy 0.0
$selall set user 0.0
$selall delete

set helix_code_list [list 1 2 3 4]

#should be defined in calling script
#set backbone_selstr "name BB" ;#martini
#set backbone_selstr "alpha" ; #charmm


#should be defined in calling script
#set protein_selstr "name BB SC1 to SC4" ;#martini
#set protein_selstr "protein" ; #charmm

#set proBeads [atomselect top $pro_selstr]
#set nProtein [expr [$proBeads num]/3892]
#$proBeads delete


set M1_start    192 
set M1_end	    216
set M2_start    218
set M2_end	243
set M3_start	251
set M3_end	272
set M4_start	286
set M4_end	307

;# Start and end residue numbers for M1 to M4
set M1_start_list [[atomselect top "resid $M1_start and $backbone_selstr"] get residue]
set M1_end_list   [[atomselect top "resid $M1_end and $backbone_selstr"] get residue]
set M2_start_list [[atomselect top "resid $M2_start and $backbone_selstr"] get residue]
set M2_end_list   [[atomselect top "resid $M2_end and $backbone_selstr"] get residue]
set M3_start_list [[atomselect top "resid $M3_start and $backbone_selstr"] get residue]
set M3_end_list   [[atomselect top "resid $M3_end and $backbone_selstr"] get residue]
set M4_start_list [[atomselect top "resid $M4_start and $backbone_selstr"] get residue]
set M4_end_list   [[atomselect top "resid $M4_end and $backbone_selstr"] get residue]

puts $M1_start_list
puts $M4_end_list
foreach start $M1_start_list end $M4_end_list chainId [list A B C D E] {
    set tempsel [atomselect top "residue $start to $end"]
    $tempsel set chain $chainId
}


[atomselect top "resid $M1_start to $M1_end and $protein_selstr"] set occupancy 1
[atomselect top "resid $M2_start to $M2_end and $protein_selstr"] set occupancy 2
[atomselect top "resid $M3_start to $M3_end and $protein_selstr"] set occupancy 3
[atomselect top "resid $M4_start to $M4_end and $protein_selstr"] set occupancy 4

