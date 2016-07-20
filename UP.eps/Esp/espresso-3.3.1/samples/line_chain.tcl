puts " "
puts "======================================================="
puts "=               My script 1: line_chain.tcl               ="
puts "======================================================="
puts " "

puts "Program Information: \n[code_info]\n"

#############################################################
#  Parameters                                               #
#############################################################

# System identification: 
set name  "line_chain"
set ident "_s1"

# System parameters
#############################################################

# Number of polymers
set n_polymers 1
# length of polymers
set l_polymers 50
# system density What is right density?! p=M/V => 50/8000
set density 0.00625

puts [format "%d %d" $n_polymers $l_polymers ]

# Interaction parameters
#############################################################

# repulsive Lennard Jonesset lj1_eps     1.0
set lj1_eps     1.0
set lj1_sig     1.0
set lj1_cut     1.12246
set lj1_shift   [calc_lj_shift $lj1_sig $lj1_cut]
set lj1_off     0.0

# attractive FENE fene_k?!
#set fene_k      30.0 
set fene_r      2.0

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
set int_n_times  5000

# More general Parameters
#############################################################

set tcl_precision 6
set vmd_output    "no"

# Initial Bond length
set bond_length   1.0

puts "Simulate the following chain in good solution:"
puts "number of polymers       : $n_polymers"
set n_monomers [expr ($n_polymers*$l_polymers)]
puts "number of monomers       : $n_monomers"
set box_length [expr pow(($n_monomers/$density),(1.0/3.0))]
puts "cubic box length         : $box_length"

# prepare observable output
set obs_file [open "$name$ident.obs" "w"]
puts $obs_file "\#$name$ident: Observables"
puts $obs_file "\#Time     R_E         R_G         E_TOT       E_KIN       E_POT"
analyze set chains 0 $n_polymers $l_polymers
set re_full [analyze re]
set re [lindex $re_full 0]
set j 0
puts "\nStart Main integration: $int_n_times times $int_steps steps"
for {set i 0} { $i <= $int_n_times } { incr i} {
    set time [setmd time]
    puts -nonewline "run $i at time=$time, R_E = $re_full (soll 15.0, deviation [expr $re/15.0])\r"
    flush stdout

    integrate $int_steps

    set f [open "$name$ident\_dist.dat" w]
    puts -nonewline $f "\#"
    puts $f "[analyze distribution 2 {0 1} 1.0 300.0 30 1 1]"
    flush $f
    close $f

    set energy [analyze energy]
    set re_full [analyze re]
    set re [lindex $re_full 0]
    set rg [lindex [analyze rg] 0]
    #set rh [lindex [analyze rh] 0]

    puts $obs_file [format "%.3e %.5e %.5e %.5e %.5e %.5e %.5e" $time $re $rg [lindex [lindex $energy 0] 1] [lindex [lindex $energy 1] 1] [lindex [lindex $energy [expr [llength $energy]-1]] 1] ]
    flush $obs_file
}
