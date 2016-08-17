puts " "
puts "======================================================="
puts "=             My script 1: chain_in_box.tcl           ="
puts "======================================================="
puts " "

puts "Program Information: \n[code_info]\n"
#  Parameters  
#############################################################
# System identification: 
set name  "chain_in_box"
set ident "_s0"
set chain_length 20

# Interaction parameters
###
set n_part_types 2

# repulsive Lennard Jonesset
#############################################################
set lj1_eps     1.0
set lj1_sig     1.0
set lj1_cut     1.12246
set lj1_shift   [calc_lj_shift $lj1_sig $lj1_cut]
set lj1_off     0.0
# attractive FENE
set fene_k      7.0 
set fene_cut    2.0

# Integration parameters
#############################################################
setmd time_step 0.0125
setmd skin      0.4
thermostat langevin 1.0 1.0
# warmup integration (with capped LJ potential)
set warm_steps   200
set warm_n_times 30
# do the warmup until the particles have at least the distance min__dist
set min_dist     0.9
# integration
set int_steps    200
set int_n_times  200

# More general Parameters
#############################################################
set tcl_precision 10
# Initial Bond length
set bond_length   1.0

puts "Simulate the chain_in_box in good solution:"
# simulation box
set box_length 7
setmd box_l $box_length $box_length $box_length 
puts "box: [setmd box_l]"

# bonded interactions
set fene_type 0;
inter $fene_type fene $fene_k $fene_cut

# non-bonded LJ interactions between all particle types
# Warning: make sure that your n_part_types really refers to all particle types that you are going to create in the next lines
for {set type1 0} { $type1 < $n_part_types } { incr type1 } {
    for {set type2 0} {  $type2 < $n_part_types } { incr type2 } {
	inter $type1 $type2 lennard-jones $lj1_eps 1.0 $lj1_cut auto
    }
}

# Particle setup
#############################################################
puts "Setup wall..."
constraint wall normal 1 1 0 type 1 penetrable 0
puts "Setup Particles (wait...)"
# set up the first particle and fix it
set i 0;
set posz 0
part 0 pos [expr $box_length/2] [expr $box_length/2] $posz  type 0 fix
puts "Created particle $i: [part $i]"
# set up the rest of the chain
set posx [expr $box_length/2]
set posy [expr $box_length/2]
set part_id 1; set mode "SAW"; set type_FENE 0  
for {set i 1} {$i < $chain_length} {incr i} {
	set theta [expr     [PI]*[t_random]]
	set phi   [expr 2.0*[PI]*[t_random]]
	set posx  [expr $posx + $bond_length*sin($theta)*cos($phi)]
	set posy  [expr $posy + $bond_length*sin($theta)*sin($phi)]
	set posz  [expr $posz + $bond_length*cos($theta)]
	part $i pos $posx $posy [expr abs($posz)] type 0 bond $fene_type [expr $i-1]
	puts "Created particle $i: [part $i]"
}

puts "part: [part]"
puts "inter: [inter]"
puts "n_part_types: [setmd n_part_types]"


set f [open "text.vtf" "w"]
writevsf $f;
writevcf $f;
close $f;

exit;


