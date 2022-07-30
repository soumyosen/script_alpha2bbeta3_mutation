
set chainA_residue_list [list 1032 771 650 596 63] 
set chainA_resname_list [list "E" "V" "T" "I" "R"]
set chainA_mutate_resname_list [list "W" "L" "M" "T" "C"]

set chainB_residue_list [list 189 433 654 662]
set chainB_resname_list [list "P" "P" "E" "R"]
set chainB_mutate_resname_list [list "S" "A" "K" "C"]

set chain_name [list "A" "B"]

set file_index 0 
############### This portion find out the interacting residues of the specified residues in chain Alpha

foreach i $chainA_residue_list res1 $chainA_resname_list res2 $chainA_mutate_resname_list {
################## for wild type        
	set close_atoms_chainA [atomselect 0 "(chain A and not resid $i) and noh and within 4.0 of (chain A and resid $i and noh)"]
        set close_atoms_chainB [atomselect 0 "chain B and noh and within 4.0 of (chain A and resid $i and noh)"]
	set resids_chainA [lsort -unique -int [$close_atoms_chainA get resid]]
        set resids_chainB [lsort -unique -int [$close_atoms_chainB get resid]]
	set close_residue_list {}
	foreach j1 $resids_chainA {
		lappend close_residue_list A_$j1
	}
        foreach j2 $resids_chainB {
                lappend close_residue_list B_$j2
        }
	puts "A$i : $close_residue_list"

#################### for mutation
        set close_atoms_chainA1 [atomselect 1 "(chain A and not resid $i) and noh and within 4.0 of (chain A and resid $i and noh)"]
        set close_atoms_chainB1 [atomselect 1 "chain B and noh and within 4.0 of (chain A and resid $i and noh)"]
        set resids_chainA1 [lsort -unique -int [$close_atoms_chainA1 get resid]]
        set resids_chainB1 [lsort -unique -int [$close_atoms_chainB1 get resid]]
        set close_residue_list1 {}
        foreach j1 $resids_chainA1 {
                lappend close_residue_list1 A_$j1
        }
        foreach j2 $resids_chainB1 {
                lappend close_residue_list1 B_$j2
        }
        puts "A$i : $close_residue_list1"
#############################################################################################################

	$close_atoms_chainA delete
        $close_atoms_chainB delete
        $close_atoms_chainA1 delete
        $close_atoms_chainB1 delete

##################### Write the output of analysis
        set new_list [lsort -unique [concat $close_residue_list $close_residue_list1]]
        puts "A$i : $new_list"

        set fp [open "contact_$file_index.dat" w]
	puts $fp "Residues,A_$res1$i,A_$res2$i"
	foreach res $new_list {
		set cond1 [lsearch $close_residue_list $res]
		set cond2 [lsearch $close_residue_list1 $res]
		if { ($cond1!=-1) && ($cond2!=-1) } {
		      puts $fp "$res,1,1"
                } elseif { ($cond1!=-1) && ($cond2==-1) } {
		      puts $fp "$res,1,0"
	        } else {
	              puts $fp "$res,0,1"
                }
        }
        flush $fp
        close $fp
        set file_index [expr $file_index+1]	

}

############### This portion find out the interacting residues of the specified residues in chain Beta

foreach i $chainB_residue_list res1 $chainB_resname_list res2 $chainB_mutate_resname_list {
############## For wild type        
	set close_atoms_chainB [atomselect 0 "(chain B and not resid $i) and noh and within 4.0 of (chain B and resid $i and noh)"]
        set close_atoms_chainA [atomselect 0 "chain A and noh and within 4.0 of (chain B and resid $i and noh)"]
	set resids_chainA [lsort -unique -int [$close_atoms_chainA get resid]]
        set resids_chainB [lsort -unique -int [$close_atoms_chainB get resid]]
	set close_residue_list {}
	foreach j1 $resids_chainA {
		lappend close_residue_list A_$j1
	}
        foreach j2 $resids_chainB {
                lappend close_residue_list B_$j2
        }
	puts "B$i : $close_residue_list"

############## For mutations
        set close_atoms_chainB1 [atomselect 1 "(chain B and not resid $i) and noh and within 4.0 of (chain B and resid $i and noh)"]
        set close_atoms_chainA1 [atomselect 1 "chain A and noh and within 4.0 of (chain B and resid $i and noh)"]
        set resids_chainA1 [lsort -unique -int [$close_atoms_chainA1 get resid]]
        set resids_chainB1 [lsort -unique -int [$close_atoms_chainB1 get resid]]
        set close_residue_list1 {}
        foreach j1 $resids_chainA1 {
                lappend close_residue_list1 A_$j1
        }
        foreach j2 $resids_chainB1 {
                lappend close_residue_list1 B_$j2
        }
        puts "B$i : $close_residue_list1"
##################################################################################################
        
        $close_atoms_chainA delete
        $close_atoms_chainB delete
        $close_atoms_chainA1 delete
        $close_atoms_chainB1 delete

##################### Write the output of analysis
        
        set new_list [lsort -unique [concat $close_residue_list $close_residue_list1]]
        puts "B$i : $new_list"

        set fp [open "contact_$file_index.dat" w]
        puts $fp "Residues,B_$res1$i,B_$res2$i"
        foreach res $new_list {
                set cond1 [lsearch $close_residue_list $res]
                set cond2 [lsearch $close_residue_list1 $res]
                if { ($cond1!=-1) && ($cond2!=-1) } {
                      puts $fp "$res,1,1"
                } elseif { ($cond1!=-1) && ($cond2==-1) } {
                      puts $fp "$res,1,0"
                } else {
                      puts $fp "$res,0,1"
                }
        }
        flush $fp
        close $fp
        set file_index [expr $file_index+1]
}
