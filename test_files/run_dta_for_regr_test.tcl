set scriptDir [file dirname [file normalize [info script]]]
set DTA_path [file join $scriptDir "../TCL/polarDensity_for_DTA.tcl"]
source $DTA_path


proc load_and_run_test {trajpath groname xtcname config home} {
	cd $trajpath
	mol new $groname
	if {$xtcname != 0} {
		mol addfile $xtcname waitfor all
		animate delete beg 0 end 0 skip 0 top
	}
	polarDensityBin $config
	mol delete top
	cd $home
}

;# Put your tests below

;# ELIC in 95% POPC 5% Cardiolipin membrane
set path [file join $scriptDir "MD_files/rep1/test_vals"]
load_and_run_test $path ../example.gro ../example.xtc ${path}/../config_for_regr_test.tcl $scriptDir

exit
