
puts " "
puts "======================================================="
puts "=       Sample script 3: one_liner_chain.tcl              ="
puts "======================================================="
puts " "

puts "Program Information: \n[code_info]\n"

#############################################################
#  Parameters                                               #
#############################################################

# System identification: 
set name  "one_liner_chain"
set ident "_s1"

# System parameters
#############################################################

# Number of polymers
set n_polymers 1
# length of polymers
set l_polymers 50
# distance between charged monomers
##set cm_distance 3
# system density
set density 0.00625

puts [format "%d %d" $n_polymers $l_polymers]

# Interaction parameters
#############################################################

# Lennard Jones
set n_part_types  1
# repulsive LJ
set lj1_cut       1.12246
set lj1_epsilon   1.0
set lj1_sig       1.0

# FENE
set fene_cut      2.0
set fene_k        7.0

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
set int_n_times  100


# More general Parameters
#############################################################

set tcl_precision 6
set vmd_output    "no"

# Initial Bond length
set bond_length   1.0

puts "number of polymers       : $n_polymers"
set n_monomers [expr ($n_polymers*$l_polymers)]
puts "number of monomers       : $n_monomers"
set npart $n_monomers
puts "total number of particles: $npart"
set box_length [expr pow(($npart/$density),(1.0/3.0))]
puts "cubic box length         : $box_length"

# Dependent Parameters
set lj1_shift  [calc_lj_shift $lj1_sig $lj1_cut]

####vtf-files#####
set vcf_file [open "$name$ident.vcf" "w"]
puts $vcf_file "\#$name$ident: coordinates"
writevcf $vcf_file
flush $vcf_file
close $vcf_file
##################


# Setup
#############################################################

# simulation box
setmd box_l $box_length $box_length $box_length 

# fene
inter 0 fene $fene_k $fene_cut

# pairwise lennard_jones for all particles
for {set ia1 0} { $ia1 < $n_part_types } { incr ia1 } {
    for {set ia2 0} { $ia2 < $n_part_types } { incr ia2 } {
	inter $ia1 $ia2 lennard-jones $lj1_epsilon 1.0 $lj1_cut $lj1_shift 0.0
    }
}

puts "Setup Particles (wait...)"
# polymers
polymer $n_polymers $l_polymers $bond_length mode PSAW  FENE 0
#puts "[part]"
set n_part [setmd n_part]
for { set i 0 } { $i < $n_part } { incr i } {
    puts "[part $i]"
}

set act_min_dist [analyze mindist]
puts "Placed [setmd n_part] particles with minimal distance: $act_min_dist"

##############
set starttime [clock seconds]


#############################################################
#  Warmup Integration                                       #
#############################################################

puts "\nStart warmup integration:"
puts "At maximum $warm_n_times times $warm_steps steps"
puts "Stop if minimal distance is larger than $min_dist"

# set LJ cap
set cap 20
inter forcecap $cap 
set i 0
while { $i < $warm_n_times && $act_min_dist < $min_dist } {
    set time [format "%8f" [setmd time]]
    puts -nonewline "run $i at time=$time (LJ cap=$cap)\r "
    flush stdout

    puts [analyze energy]
    integrate $warm_steps

    if { $vmd_output=="yes" } { imd positions }

    set act_min_dist [analyze mindist]
    puts -nonewline "minimal distance = $act_min_dist\r"
    flush stdout
#   Increase LJ cap
    set cap [expr $cap+10]
    inter forcecap $cap
    incr i
}

puts "\nWarmup done: Minimal distance [analyze mindist]"
puts "Verlet reuses: [setmd verlet_reuse]"

#############################################################
#      Integration                                          #
#############################################################

#
#      Write blockfiles for restart
#############################################################
set trajectory [open "$name$ident.config" "w"]

blockfile $trajectory write variable {box_l time_step skin}
blockfile $trajectory write interactions
blockfile $trajectory write integrate
blockfile $trajectory write thermostat
flush $trajectory

# prepare observable output
set obs_file [open "$name$ident.obs" "w"]
puts $obs_file "\#$name$ident: Observables"
puts $obs_file "\#Time     Cputime     R_E         R_G          E_TOT       E_KIN       E_P"
analyze set chains 0 $n_polymers $l_polymers
set re_full [analyze re]
set re [lindex $re_full 0]
set j 0
puts "\nStart Main integration: $int_n_times times $int_steps steps"
for {set i 0} { $i <= $int_n_times } { incr i} {
    set time [setmd time]
    #~ set timingstart [clock clicks -milliseconds]
    puts -nonewline "run $i at time=$time, R_E = $re_full (soll 15.0, deviation [expr $re/15.0])\r"
    flush stdout

    integrate $int_steps
    set timingcurr [clock clicks -milliseconds]
    
    
    set elapsedtime [expr  $timingcurr - $timingstart]
    #puts  "elapsed time: $elapsedtime"

    set f [open "$name$ident\_dist.dat" w]
    puts -nonewline $f "\#"
    puts $f "[analyze distribution 2 {0 1} 1.0 300.0 30 1 1]"
    flush $f
    close $f

    set energy [analyze energy]
    set re_full [analyze re]
    set re [lindex $re_full 0]
    set rg [lindex [analyze rg] 0]
 
    set vsf_file [open "$name$ident.vsf" "w"] 
    puts $vsf_file "\#$name$ident: topology"
    writevsf $vsf_file
    flush $vsf_file
    close $vsf_file   
 
 puts $obs_file [format "%.3e %.3e %.5e %.5e %.5e %.5e %.5e" $time $elapsedtime $re $rg [lindex [lindex $energy 0] 1] [lindex [lindex $energy 1] 1] [expr [analyze energy total] - [analyze energy kinetic]] ]
    flush $obs_file
    if { $vmd_output=="yes" } { imd positions }
    #   write intermediate configuration
    if { $i%50==0 } {
	blockfile $trajectory write particles
	flush $trajectory
	incr j
    }
}

puts "\nIntegration done."

#############################################################

puts "\nFinished"
exit
