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
############# for protein indexes and another list of water indexes. 
proc index_in_hbonds {molid wat_residue_list ch resid} {
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
		        	} else {
		        	        lappend water_index $ind1
		        	}
		        }	
		} else {
		        foreach ind1 $donors_or_acceptors1 {
		        	if {$ind1<=$index_cut1} {
		        		lappend protein_index $ind1
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
	set index_lists [list $protein_index $water_index]
	return $index_lists
}   

############## Take the protein indexes which make Hbonds with the first layer of water and return the information of those
############## protein residues. It also return another list which used later to remove the repetation of the same residues
############## in two water mediated interactions
proc single_wat_mediated {molid list_prot_index} {
        set all_residue_info {}
        set repeat_info {}	
	foreach ind $list_prot_index {
		set sel_residue [atomselect $molid "index $ind"]
		set residue_info {}
		set residue_chain [$sel_residue get chain]
		set residue_name [$sel_residue get resname]
		set residue_num [$sel_residue get resid]
		set comb_chain_num $residue_chain$residue_num
		lappend residue_info "$residue_chain $residue_name $residue_num 1water"
		if {[lsearch $repeat_info $comb_chain_num] == -1} {
			lappend all_residue_info [lindex $residue_info 0]
			lappend repeat_info $comb_chain_num 
		}
		$sel_residue delete
	}
	set residue_info_list [list $all_residue_info $repeat_info]
	return $residue_info_list
}	

############## Take the protein indexes which make Hbonds with the 2nd layer of water and return the information of the protein 
############## residues. The information is appended into the list of 1 water mediated interacting protein residues. If a
############## residue already makes interaction with the desired residue by single water mediated that residue will not be included
############## in this list
proc two_wats_mediated {molid list_prot_index all_residue_info repeat_info} {
        foreach ind $list_prot_index {
                set sel_residue [atomselect $molid "index $ind"]
                set residue_info {}
                set residue_chain [$sel_residue get chain]
                set residue_name [$sel_residue get resname]
                set residue_num [$sel_residue get resid]
		set comb_chain_num $residue_chain$residue_num
                lappend residue_info "$residue_chain $residue_name $residue_num 2waters"
                if {[lsearch $repeat_info $comb_chain_num] == -1} {
			lappend all_residue_info [lindex $residue_info 0]
                        lappend repeat_info $comb_chain_num
		}
                $sel_residue delete
        }
        
	set residue_info_list [list $all_residue_info $repeat_info]
	#set residue_info_list $all_residue_info
        return $residue_info_list
}

############# Main work is done in this function which used the previous functions. Take the information
############# of the chain and resid number of the desired residues 
proc results {chain resid} {

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

	set d_or_a_wat_index0 [select_water_indexes 0 $donors_or_acceptors0]
	set d_or_a_wat_index1 [select_water_indexes 1 $donors_or_acceptors1]
	set d_or_a_wat_residues0 [find_water_residues 0 $d_or_a_wat_index0]
	set d_or_a_wat_residues1 [find_water_residues 1 $d_or_a_wat_index1]

        set all_indexes_1water0 [index_in_hbonds 0 $d_or_a_wat_residues0 $chain $resid]
        set all_indexes_1water1 [index_in_hbonds 1 $d_or_a_wat_residues1 $chain $resid]
        
	set single_water_int0 [single_wat_mediated 0 [lindex $all_indexes_1water0 0]]
        set single_water_int1 [single_wat_mediated 1 [lindex $all_indexes_1water1 0]]

	set all_1water_res0 [find_water_residues 0 [lindex $all_indexes_1water0 1]]
	set all_1water_res1 [find_water_residues 1 [lindex $all_indexes_1water1 1]]
        
        set all_indexes_2waters0 [index_in_hbonds 0 $all_1water_res0 $chain $resid]
        set all_indexes_2waters1 [index_in_hbonds 1 $all_1water_res1 $chain $resid]
         
        set interaction_info0 [lindex $single_water_int0 0]
        set checking_info0 [lindex $single_water_int0 1]
        set interaction_info1 [lindex $single_water_int1 0]
        set checking_info1 [lindex $single_water_int1 1]
        
	#puts "single_water_int0 $single_water_int0"
	#puts "molid 0: protein index [lindex $all_indexes_2waters0 0] full previous information $interaction_info0 check info $checking_info0"
	#puts "single_water_int1 $single_water_int1"
	#puts "molid 1: protein index [lindex $all_indexes_2waters1 0] full previous information $interaction_info1 check info $checking_info1"
        set two_waters_int0 [two_wats_mediated 0 [lindex $all_indexes_2waters0 0] $interaction_info0 $checking_info0 ]
        set two_waters_int1 [two_wats_mediated 1 [lindex $all_indexes_2waters1 0] $interaction_info1 $checking_info1 ]
	#puts "two_waters_int0 $two_waters_int0"
	#puts "two_waters_int1 $two_waters_int1"
        #puts "molid 0 $two_waters_int0"
	#puts "molid 1 $two_waters_int1"
	set results_list [list $two_waters_int0 $two_waters_int1]
        #puts $results_list
        return $results_list 
}

proc write_output {file_index interaction_info int_info_wild int_info_mutate wild_type mutation_type} {
	set fil [open "water_mediate_$file_index.dat" w]
	puts $fil "protein,chain,res_num,amino_acid,int_type,$wild_type,$mutation_type"
        if {[llength $interaction_info]!=0} {
        	foreach interaction $interaction_info {
        		set chain [lindex $interaction 0]
        	        set res_name [lindex $interaction 1]
        	        set res_num [lindex $interaction 2]
        	        set int_type [lindex $interaction 3]
        	        set protein $res_name$res_num.$chain	      
                        set cond0 [lsearch $int_info_wild $interaction]
                        set cond1 [lsearch $int_info_mutate $interaction]
        		if { ($cond0!=-1) && ($cond1!=-1) } { 
        		        set show0 True
        		        set show1 True
        		} elseif { ($cond0!=-1) && ($cond1==-1) } {
        		        set show0 True
        		        set show1 False
        	        } else {
        		        set show0 False
        		        set show1 True
        		}
                        puts $fil "$protein,$chain,$res_num,$res_name,$int_type,$show0,$show1"
                        #puts "$protein,$chain,$res_num,$res_name,$int_type,$show0,$show1"
        		}
        
        	}
	flush $fil
	close $fil
}

proc append_info_rep {list_where_append appending_list0 appending_list1} {
	foreach elem0 $appending_list0 {
		lappend list_where_append $elem0
	}
	foreach elem1 $appending_list1 {
		lappend list_where_append $elem1
	}
	return $list_where_append
}



############# Final execution
set file_index 0
foreach resid $chainA_residue_list {
	puts "A $resid"
	set res_name0 [lsort -unique [[atomselect 0 "chain A and resid $resid and sidechain"] get resname]]
	set res_name1 [lsort -unique [[atomselect 1 "chain A and resid $resid and sidechain"] get resname]]
	set info0 $res_name0$resid.A
	set info1 $res_name1$resid.A
	set result [results A $resid]
	puts "full result $result"
	set int_info0 [lindex [lindex $result 0] 0]
	puts "int_info0 $int_info0"
	set int_rep0 [lindex [lindex $result 0] 1]
	puts "int_rep0 $int_rep0"
	set int_info1 [lindex [lindex $result 1] 0]
	puts "int_info1 $int_info1"
	set int_rep1 [lindex [lindex $result 1] 1]  
	puts "int_rep1 $int_rep1"
	#puts "A $resid $result"
	set int_info01 {}
	set int_rep01 {}
        set int_info01_a [lsort -unique [append_info_rep $int_info01 $int_info0 $int_info1]]
	set int_rep01_a [lsort -unique [append_info_rep $int_rep01 $int_rep0 $int_rep1]]
	puts "int_info01_a $int_info01_a"
	puts "int_rep01_a $int_rep01_a"
	write_output $file_index $int_info01_a $int_info0 $int_info1 $info0 $info1
	set file_index [expr $file_index+1]
}


foreach resid $chainB_residue_list {
	puts "B $resid"
	set res_name0 [lsort -unique [[atomselect 0 "chain B and resid $resid and sidechain"] get resname]]
	set res_name1 [lsort -unique [[atomselect 1 "chain B and resid $resid and sidechain"] get resname]]
	set info0 $res_name0$resid.B
	set info1 $res_name1$resid.B
	set result [results B $resid]
	puts "full result $result"
	set int_info0 [lindex [lindex $result 0] 0]
	puts "int_info0 $int_info0"
	set int_rep0 [lindex [lindex $result 0] 1]
	puts "int_rep0 $int_rep0"
	set int_info1 [lindex [lindex $result 1] 0]
	puts "int_info1 $int_info1"
	set int_rep1 [lindex [lindex $result 1] 1]  
	puts "int_rep1 $int_rep1"
	#puts "B $resid $result"
	set int_info01 {}
	set int_rep01 {}
        set int_info01_a [lsort -unique [append_info_rep $int_info01 $int_info0 $int_info1]]
	set int_rep01_a [lsort -unique [append_info_rep $int_rep01 $int_rep0 $int_rep1]]
	puts "int_info01_a $int_info01_a"
	puts "int_rep01_a $int_rep01_a"
	write_output $file_index $int_info01_a $int_info0 $int_info1 $info0 $info1
	set file_index [expr $file_index+1]
}


