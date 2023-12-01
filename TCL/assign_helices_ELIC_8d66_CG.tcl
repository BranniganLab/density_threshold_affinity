set selall [atomselect top {all}]
$selall set occupancy 0.0
$selall set user 0.0

set M1start    202 
set M1end	226
set M2start    228
set M2end	253
set M3start	261
set M3end	282
set M4start	296
set M4end	317

set M1 [atomselect top "resid $M1start to $M1end and name BB SC1 to SC4"]
$M1 set occupancy 1
$M1 delete

set M2 [atomselect top "resid $M2start to $M2end and name BB SC1 to SC4"] 
$M2 set occupancy 2
$M2 delete

set M3 [atomselect top "resid $M3start to $M3end and name BB SC1 to SC4"] 
$M3 set occupancy 3
$M3 delete

set M4 [atomselect top "resid $M4start to $M4end and name BB SC1 to SC4"] 
$M4 set occupancy 4
$M4 delete

