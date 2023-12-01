
;# I found out the acyl-chains for the PAP1 to 3 are backwards, this should fix it
set test [atomselect top "resname PAP1 PAP2 PAP3"]
if {[$test num] > 0} {
	$test delete
	set org [list C1A C2A C3A C4A D1B D2B D3B D4B C5B]
	set chg [list C1B C2B C3B C4B D1A D2A D3A D4A C5A]
	foreach o $org c $chg {
		set sel [atomselect top "resname PAP1 PAP2 PAP3 and name $o"]
		$sel set name $c
		$sel delete
	}
} else {$test delete}
# TODO These need a "Holding" variable to prevent being overwritten
set test [atomselect top "resname POSM"]
if {[$test num] > 0} {
	$test delete
	set org [list T1A C2A C3A C1B D2B C3B C4B T1C C2C C3C]
	set chg [list T1C C2C C3C C1A D2A C3A C4A T1B C2B C3B]
	foreach o $org c $chg {
		set sel [atomselect top "resname POSM and name $o"]
		$sel set name $c
		$sel delete
	}
} else {$test delete}

set test [atomselect top "resname PNSM"]
if {[$test num] > 0} {
	$test delete
	set org [list T1A C2A C3A C1B C2B C3B D4B C5B C6B T1C C2C C3C]
	set chg [list T1C C2C C3C C1A C2A C3A D4A C5A C6A T1B C2B C3B]
	foreach o $org c $chg {
		set sel [atomselect top "resname PNSM and name $o"]
		$sel set name $c
		$sel delete
	}
} else {$test delete}

;#Lipid Stuff
atomselect macro lipids {all and (not resname W ION) and (not name BB SC1 to SC4 NH3 NC3 GL1 GL2 AM1 AM2 CNO CN0 GL0 GLO) and not (name C1 C2 C3 P1 P2 P3  and resname ".*PI" ".*P1" ".*P2" ".*P3")}
atomselect macro PC {lipids and resname ".*PC"}
atomselect macro PE {lipids and resname ".*PE"}
atomselect macro PS {lipids and resname ".*PS"}
atomselect macro PA {lipids and resname ".*PA"}
atomselect macro PI {lipids and resname ".*PI"}
atomselect macro SM {lipids and resname ".*SM"}
atomselect macro PG {lipids and resname ".*PG"}
atomselect macro P1a {lipids and resname ".*P1"}
atomselect macro P2a {lipids and resname ".*P2"}
atomselect macro P3a {lipids and resname ".*P3"}

;#Sat
atomselect macro din0 {lipids and (resname "DL.*" "DP.*" "DX.*")  and (not name AM1 AM2)}
atomselect macro sn0 {lipids and resname "P.*" "X.*" and (name ".*B" PO4) and not name BB} 
atomselect macro n0 {lipids din0 or sn0}

;#Mono
atomselect macro din9 {lipids and (resname "DO.*" "DV.*" "DG.*" "OG.*" "ON.*") }
atomselect macro sn9 {resname "PO.*" "PV.*" "PG.*" "PN.*" and name ".*A" PO4  and not name BB} 
atomselect macro sn9b {resname "OU.*" "OI.*" "OA.*" and name ".*B" PO4 and not name BB} 

atomselect macro n9 {lipids din9 or sn9 or sn9b}

;#Di removed. These are really n-6
# atomselect macro din7 {lipids and (resname  "DF.*") }
# atomselect macro sn7 {(resname  "PF.*" "OF.*" "IF.*") and name ".*A" PO4 and not name BB} 
# atomselect macro n7 {lipids din7 or sn7}

;#n6
atomselect macro din6 {lipids and (resname "DF.*" "DE.*" "DQ.*" "DI.*" "IQ.*" "IE.*" "DA.*" "IA.*" "IF.*") }
atomselect macro sn6 {(resname "OF.*" "PF.*" "PE.*" "PQ.*" "PI.*" "OI.*" "PA.*" "OA.*") and name ".*A" PO4 and not name BB}
atomselect macro sn6b {(resname "IA.*" "IE.*" "IU.*" "IF.*") and name ".*B" PO4 and not name BB}
;#atomselect macro sn6b {(resname "IE.*") and name ".*A" PO4 and not name BB} 

atomselect macro n6 {lipids din6 or sn6 or sn6b}

;#n3
atomselect macro din3 {lipids and (resname  "DU.*") }
atomselect macro sn3 {(resname "PU.*" "IU.*" "OU.*" ) and name ".*A" PO4 and not name BB} 
atomselect macro n3 {lipids din3 or sn3}

;#chol1
atomselect macro chol {resname CHOL}

color change rgb 5 0.960000 0.880000 0.730000

;# Protien Stuff
if { [file exists ~/lms464/github/JPC_Special/common/grace/assign_helices_2BG9_CG_lms2.tcl]  } {
    set bb [atomselect top "name BB"]
    if {[$bb num] > 0} {
        catch {
            source ~/lms464/github/JPC_Special/common/grace/assign_helices_2BG9_CG_lms2.tcl
            color Chain A green3
            color Chain B purple
            color Chain C silver
            color Chain D green3
            color Chain E cyan3
        }
	}
	$bb delete
}

proc Make_Sel {sel col} {
    mol color ColorID $col
    mol representation VDW 1.5 32.0
    mol selection  $sel
    mol material AOChalky
    mol addrep top
}

proc Chains {} {

    mol delrep 0 top
    
    Make_Sel "n0" 23
    Make_Sel "n9" 3
    Make_Sel "n7" 28
    Make_Sel "n6" 13
    Make_Sel "n3" 5
    Make_Sel "resname CHOL" 1

}
#TODO set up macro for head groups


proc Head {} {
    mol delrep 0 top
    Make_Sel "PC" 23
    Make_Sel "PE" 8
    Make_Sel "PS" 28
    Make_Sel "PI" 32
    Make_Sel "PA" 5
    Make_Sel "SM" 19
    Make_Sel "PG" 14
    Make_Sel "resname CHOL" 1
}


