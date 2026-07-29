set scriptDir [file dirname [file normalize [info script]]]
set DTA_path [file normalize [file join $scriptDir "../../TCL/polarDensity_for_DTA.tcl"]]
source $DTA_path

proc load_traj {trajpath groname xtcname} {
	mol new $trajpath$groname
	if {$xtcname != 0} {
		mol addfile $trajpath$xtcname waitfor all
		animate delete beg 0 end 0 skip 0 top
	}
}

load_traj "../MD_files/rep1/" "example.gro" "example.xtc"

proc RunBenchmark {} {
	polarDensityBin "../MD_files/rep1/config_for_regr_test.tcl"
}