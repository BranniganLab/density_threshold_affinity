set script_path [pwd]
source $script_path/../TCL/polarDensity_for_DTA.tcl


proc load_and_run_test {trajpath groname xtcname config} {
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

;# E protein cartesian
set path "${script_path}/E-protein_trajectory"
load_and_run_test $path DT_test.gro DT_test.xtc ${path}/config_for_regr_test.tcl

exit
