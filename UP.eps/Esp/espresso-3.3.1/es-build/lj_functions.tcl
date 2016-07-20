#source lj_functions.tcl 
puts ""
puts "========================"
puts "=lj_liquid_tutorial.tcl="
puts "========================"
puts ""
puts "Espresso Code Base : \n[code_info]\n"
#cellsystem domain_decomposition -n_verlet _list 
# System identification:
set name "lj_liquid"
set ident "_s1"

#System parameters 
##################
## 

set n_part 108

#Interaction parameters 
#######################
##
#
set lj1_eps 1.0
set lj1-sig 1.0
set lj1_cut 2.5
set lj1_shift [expr -(pow(1.0/$lj1_cut,12)-pow(1.0/$lj1_cut,6))]
set lj1_offset 0.0

#Integration parameters
######################
#

thermostat off; 

setmd time_step 0.001
set eq_tstep 0.0001
set tstep 0.001
set skin 0.1
setmd skin $skin
set target_temperature 0.728

##

set warm_steps 100
set warm_n_times 2000

#
#

set min_dist 0.87

#

set sampling_interval 1000
set equilibration_interval 1000

set sampling_interations 200
set equilibration_iterations 200


#Other parameters
##

set tcl_precision 8
#
expr srand([pid])
#Particle setup
#

set density 0.8442

#
set box_length [expr pow ($n_part/$density,1.0/3.0)+2*$skin]
puts "density = $density box_length = $box_length"
setmd box $box_length $box_length $box_length

#
for { set i 0} { $i < $n_part } { incr i} {
	set pos_x [expr rand()*$box_length]
	set pos_y [expr rand()*$box_length]
	set pos_z [expr rand()*$box_length]
	part $i pos $pos_x $pos_y $pos_z q 0.0 type 0
}

#
writepdb data/config.pdb

#Interaction setup
##############
##
##

inter 0 0 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift $lj1_offset
