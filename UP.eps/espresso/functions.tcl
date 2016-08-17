### This function will be set up posx/posy.
proc position_maker {box_lenght} {
puts "Not yet."
}

### This function will be produce polymer chain.
proc make_polymer {my_index chain_length bond_length fene_type posx posy} {
#set up the first particle and fix it 
#set i $my_index
set posz 0
part $my_index pos $posx $posy $posz type 0 fix
puts "Created particle $my_index: [part $my_index]"
# set up the rest of the chain
#set part_id 1; set mode "SAW"; set type_FENE 0
for {set i [expr $my_index+1]} {$i < [expr $chain_length+$my_index]} {incr i} {
        set theta [expr     [PI]*[t_random]]
        set phi   [expr 2.0*[PI]*[t_random]]
        set posx  [expr $posx + $bond_length*sin($theta)*cos($phi)]
        set posy  [expr $posy + $bond_length*sin($theta)*sin($phi)]
        set posz  [expr $posz + [expr abs($bond_length*cos($theta))]]
	part $i pos $posx $posy $posz type 0 bond $fene_type [expr $i-1]
	puts "Created particle $i: [part $i]"

	#set f [open "text.vtf" "w"]
        #writevsf $f;
        #writevcf $f;
        #close $f;
}
}

