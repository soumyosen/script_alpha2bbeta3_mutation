set AA_3letters_code [list "ARG" "HIS" "LYS" "ASP" "GLU" "SER" "THR" "ASN" "GLN" "CYS" "SEC" "GLY" "PRO" "ALA" "VAL" "ILE" "LEU" "MET" "PHE" "TYR" "TRP"]
set AA_1letter_code [list "R" "H" "K" "D" "E" "S" "T" "N" "Q" "C" "U" "G" "P" "A" "V" "I" "L" "M" "F" "Y" "W"]

set chainA_residue_list [list 771 650 596 456]
set chainA_resname_list [list "V" "T" "I" "V"]
set chainA_mutate_resname_list [list "L" "M" "T" "I"]

set chainB_residue_list [list 63 119 189 228 329 662]
set chainB_resname_list [list "R" "R" "P" "R" "N" "R"]
set chainB_mutate_resname_list [list "C" "Q" "S" "H" "D" "C"]

set index_cut0 26863
set index_cut1 26822
########### Function to find out the indexes of water. 26863 is the last index of wild type protein system. 26846 is the last index of mutatedprotein system
proc select_water_indexes {molid list_of_index} {
        set water_indexes {}
        global index_cut0 index_cut1 
	if {$molid==0} {
	        foreach ind $list_of_index {
		        if $ind>$index_cut0 {
		                lappend water_indexes $ind
	                } 
                }
	} else {
		foreach ind $list_of_index {
                        if $ind>$index_cut1 {
                                lappend water_indexes $ind
                        }
                }
	}
        return $water_indexes
}

########### Function to get the water residue numbers using the water indexes. Take out only the unique residues
proc find_water_residues {molid list_of_water_index} {
	set water_residues {}
        foreach ind $list_of_water_index {
		set sel1 [atomselect $molid "index $ind"]
	        set res [$sel1 get residue]
	        lappend water_residues $res
	        $sel1 delete
        }
        set water_residues [lsort -unique $water_residues]
        return $water_residues
}

############# Find out the list of indexes which are in hbonds with a set of water residues. Return a list of lists. One list is 
############# for protein indexes and another list of water indexes and another list of water residues which are making Hbonds with protein. 
proc index_in_hbonds {molid wat_residue_list ch resid water_res} {
	set protein_index {}
	set water_index {}
        global index_cut0 index_cut1 
	foreach res $wat_residue_list {
		set sel4 [atomselect $molid "residue $res"]
		set sel5 [atomselect $molid "(not residue $res and not (chain $ch and resid $resid)) and ((water or (protein and sidechain)) and name \"N.*\" \"O.*\" \"S.*\" FA F1 F2 F3) and within 8 of (residue $res)"]
                set hbonds3 [measure hbonds 3.5 50 $sel4 $sel5]
                set hbonds4 [measure hbonds 3.5 50 $sel5 $sel4]
                set donors1 [concat [lindex $hbonds3 0] [lindex $hbonds4 0]]
		set acceptors1 [concat [lindex $hbonds3 1] [lindex $hbonds4 1]]
                set donors_or_acceptors1 [lsort -unique [concat $donors1 $acceptors1]]
                if {$molid==0} {
		        foreach ind1 $donors_or_acceptors1 {
		        	if {$ind1<=$index_cut0} {
		        		lappend protein_index $ind1
					lappend water_res $res
		        	} else {
		        	        lappend water_index $ind1
		        	}
		        }	
		} else {
		        foreach ind1 $donors_or_acceptors1 {
		        	if {$ind1<=$index_cut1} {
		        		lappend protein_index $ind1
					lappend water_res $res
		        	} else {
		        	        lappend water_index $ind1
		        	}
		        }	
		}

                $sel4 delete
		$sel5 delete
	}       
        set protein_index [lsort -unique $protein_index]
	set water_index [lsort -unique $water_index]
	set index_lists [list $protein_index $water_index $water_res]
	return $index_lists
}   


######### Finding out just interacting water residues
proc search_int_wat_residues {first_layer list_of_wat_res molid} {
	set interacting_wat_res {}
	foreach res $list_of_wat_res {
		if {[lsearch $first_layer $res]!=-1} {
			########## Water residue from the first layer
			lappend interacting_wat_res $res
		} else {
			######### Water residue from the second layer, but need to know which water residue of the first layer interacts with it.
			lappend interacting_wat_res $res
			set sel6 [atomselect $molid "residue $res"]
		        set sel7 [atomselect $molid "residue $first_layer"]
			set hbonds5 [measure hbonds 3.5 50 $sel6 $sel7]
                        set hbonds6 [measure hbonds 3.5 50 $sel7 $sel6]
                        set donors2 [concat [lindex $hbonds5 0] [lindex $hbonds6 0]]
                        set acceptors2 [concat [lindex $hbonds5 1] [lindex $hbonds6 1]]
                        set donors_or_acceptors2 [lsort -unique [concat $donors2 $acceptors2]]
	                set residues [find_water_residues $molid $donors_or_acceptors2]
                        foreach res1 $residues {
				lappend interacting_wat_res $res1
			}
			$sel6 delete
			$sel7 delete
		}
	}
	set interacting_wat_res [lsort -unique $interacting_wat_res]
	return $interacting_wat_res
}


############# Main work is done in this function which used the previous functions. Take the information
############# of the chain and resid number of the desired residues 
proc results {chain resid} {
        puts "chain $chain and resid $resid"
	set sel0_1 [atomselect 0 "(chain $chain and resid $resid and sidechain) and (name \"N.*\" \"O.*\" \"S.*\" FA F1 F2 F3)"]
        set sel0_2 [atomselect 0 "same residue as water and within 8 of (chain $chain and resid $resid and sidechain)"]
        set sel1_1 [atomselect 1 "(chain $chain and resid $resid and sidechain) and (name \"N.*\" \"O.*\" \"S.*\" FA F1 F2 F3)"]
        set sel1_2 [atomselect 1 "same residue as water and within 8 of (chain $chain and resid $resid and sidechain)"]
         
        set hbonds0_1 [measure hbonds 3.5 50 $sel0_1 $sel0_2]
        set hbonds0_2 [measure hbonds 3.5 50 $sel0_2 $sel0_1]
        set hbonds1_1 [measure hbonds 3.5 50 $sel1_1 $sel1_2]
        set hbonds1_2 [measure hbonds 3.5 50 $sel1_2 $sel1_1]
        #puts "hbonds0_1 $hbonds0_1"
        #puts "hbonds0_2 $hbonds0_2"
        #puts "hbonds1_1 $hbonds1_1"
        #puts "hbonds1_2 $hbonds1_2"

        set donors0 [concat [lindex $hbonds0_1 0] [lindex $hbonds0_2 0]]
        set acceptors0 [concat [lindex $hbonds0_1 1] [lindex $hbonds0_2 1]]
        set donors1 [concat [lindex $hbonds1_1 0] [lindex $hbonds1_2 0]]
        set acceptors1 [concat [lindex $hbonds1_1 1] [lindex $hbonds1_2 1]]
        #puts "donors0 $donors0"
        #puts "acceptors0 $acceptors0"
        #puts "donors1 $donors1"
        #puts "acceptors1 $acceptors1"


        set donors_or_acceptors0 [lsort -unique [concat $donors0 $acceptors0]]
        set donors_or_acceptors1 [lsort -unique [concat $donors1 $acceptors1]]
        #puts "donors_or_acceptors0 $donors_or_acceptors0"
        #puts "donors_or_acceptors1 $donors_or_acceptors1"


        $sel0_1 delete
        $sel0_2 delete
        $sel1_1 delete
        $sel1_2 delete
        
	set water_res0 {}
	set water_res1 {}

	set d_or_a_wat_index0 [select_water_indexes 0 $donors_or_acceptors0]
	set d_or_a_wat_index1 [select_water_indexes 1 $donors_or_acceptors1]
	
	set d_or_a_wat_residues0 [find_water_residues 0 $d_or_a_wat_index0]
	set d_or_a_wat_residues1 [find_water_residues 1 $d_or_a_wat_index1]

        set all_indexes_1water0 [index_in_hbonds 0 $d_or_a_wat_residues0 $chain $resid $water_res0]
        set all_indexes_1water1 [index_in_hbonds 1 $d_or_a_wat_residues1 $chain $resid $water_res1]
        #puts "all_indexes_1water0 $all_indexes_1water0"
        #puts "all_indexes_1water1 $all_indexes_1water1"
        #puts "lindex 0 [lindex $all_indexes_1water0 0]"
        #puts "lindex 1 [lindex $all_indexes_1water0 1]"
        #puts "lindex 2 [lindex $all_indexes_1water0 2]"


	set all_1water_res0 [find_water_residues 0 [lindex $all_indexes_1water0 1]]
	set all_1water_res1 [find_water_residues 1 [lindex $all_indexes_1water1 1]]
        
        #set all_indexes_2waters0 [index_in_hbonds 0 $all_1water_res0 $chain $resid $water_res0]
	set all_indexes_2waters0 [index_in_hbonds 0 $all_1water_res0 $chain $resid [lindex $all_indexes_1water0 2]]
        #set all_indexes_2waters1 [index_in_hbonds 1 $all_1water_res1 $chain $resid $water_res1]
	set all_indexes_2waters1 [index_in_hbonds 1 $all_1water_res1 $chain $resid [lindex $all_indexes_1water1 2]]

        #puts "all_indexes_2waters0 $all_indexes_2waters0"
        #puts "all_indexes_2waters1 $all_indexes_2waters1"
        #puts "lindex 0 [lindex $all_indexes_2waters0 0]"
        #puts "lindex 1 [lindex $all_indexes_2waters0 1]"
        #puts "lindex 2 [lindex $all_indexes_2waters0 2]"
	
	#puts "molid 0 1st layer of water around prot_residues $d_or_a_wat_residues0"
        #puts "molid 0 water residues interacting with protein 1 and 2 water mediate [lindex $all_indexes_2waters0 2]"
	set water_residues_0 [search_int_wat_residues $d_or_a_wat_residues0 [lindex $all_indexes_2waters0 2] 0]
	#puts "molid 0 final interacting residues $water_residues_0"
	#puts "molid 1 1st layer of water around prot_residues $d_or_a_wat_residues1"
        #puts "molid 1 water residues interacting with protein 1 and 2 water mediate [lindex $all_indexes_2waters1 2]"
        set water_residues_1 [search_int_wat_residues $d_or_a_wat_residues1 [lindex $all_indexes_2waters1 2] 1]
        #puts "molid 1 final interacting residues $water_residues_1"
        set results_list [list $water_residues_0 $water_residues_1]
        #puts $results_list
        return $results_list
}

proc append_water_residues {list_of_wat_residues_append wat_residues_for_resid} {
	foreach elem0 $wat_residues_for_resid {
                lappend list_of_wat_residues_append $elem0
        }
	return $list_of_wat_residues_append
}


##### final execution
set water_residue_list_0 {}
set water_residue_list_1 {}

foreach resid $chainA_residue_list {
	puts "A $resid"
	set result [results A $resid]
        puts "full result $result"
        set wat_res_0 [lindex $result 0]
	#puts "wat_res_0 $wat_res_0"
        set wat_res_1 [lindex $result 1]
	#puts "wat_res_1 $wat_res_1"
        set water_residue_list_0 [append_water_residues $water_residue_list_0 $wat_res_0]
        set water_residue_list_1 [append_water_residues $water_residue_list_1 $wat_res_1]

}

foreach resid $chainB_residue_list {
	puts "B $resid"
	set result [results B $resid]
        puts "full result $result"
        set wat_res_0 [lindex $result 0]
	#puts "wat_res_0 $wat_res_0"
        set wat_res_1 [lindex $result 1]
	#puts "wat_res_1 $wat_res_1"
        set water_residue_list_0 [append_water_residues $water_residue_list_0 $wat_res_0]
        set water_residue_list_1 [append_water_residues $water_residue_list_1 $wat_res_1]
}

puts "molid 0 final water residue list $water_residue_list_0"
puts "molid 1 final water residue list $water_residue_list_1"

set fil0 [open "wat_residues_wild1.dat" w]
puts $fil0 "$water_residue_list_0"
flush $fil0
close $fil0

set fil1 [open "wat_residues_mutate.dat" w]
puts $fil1 "$water_residue_list_1"
flush $fil1
close $fil1





