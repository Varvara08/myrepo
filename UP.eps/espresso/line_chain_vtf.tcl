puts " "
puts "======================================================="
puts "=               My script 1: line_chain.tcl           ="
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
set l_polymers 30
# system density What is right density?! p=M/V => 50/8000; box_length 20.
set density [expr 20./8000]

puts [format "%d %d" $n_polymers $l_polymers ]

# Interaction parameters
#############################################################

# Lennard Jones
set n_part_types  1
# repulsive Lennard Jonesset lj1_eps     1.0
set lj1_eps     1.0
set lj1_sig     1.0
set lj1_cut     1.12246
set lj1_shift   [calc_lj_shift $lj1_sig $lj1_cut]
set lj1_off     0.0

# attractive FENE fene_k?!
set fene_k      7.0 
set fene_cut    2.0

# Integration parameters
#############################################################
#time_step 0.0125
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
set int_n_times  500

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
set n_part $n_monomers
puts "total number of particles: $n_part"
#What is this num?! 1/3?
set box_length [expr pow(($n_part/$density),(1.0/3.0))] 
puts "box length         : $box_length"
puts "density of polymer : $density"

# Dependent Parameters
set lj1_shift  [calc_lj_shift 1 $lj1_cut]


# Creating vsf topology
set vsf_file [open ./data/m$l_polymers/$name$ident\.vsf "w"]
puts $vsf_file "\#Some info"
writevsf $vsf_file 
flush $vsf_file
close $vsf_file

# Setup
#############################################################

# simulation box
setmd box_l $box_length $box_length $box_length 

# fene
inter 0 fene $fene_k $fene_cut

# pairwise lennard_jones for all particles
for {set ia1 0} { $ia1 < $n_part_types } { incr ia1 } {
    for {set ia2 0} { $ia2 < $n_part_types } { incr ia2 } {
	inter $ia1 $ia2 lennard-jones $lj1_eps 1.0 $lj1_cut $lj1_shift 0.0
    }
}

puts "Setup Particles (wait...)"
# polymers
polymer $n_polymers $l_polymers $bond_length mode PSAW FENE 0
#puts "[part]"
set n_part [setmd n_part]
for { set i 0 } { $i < $n_part } { incr i } {
    puts "[part $i]"
}

set act_min_dist [analyze mindist]
puts "Placed [setmd n_part] particles with minimal distance: $act_min_dist"

#I think this is the need interaction
#inter 0 0 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift $lj1_off
#inter 0   FENE          $fene_k  $fene_cut

##############################################################
#  prepare vmd connection with VTF                          #
#############################################################
#    if { $vmd_output=="yes" } {
#	puts -nonewline "\nWrite vsf and vcf for VMD connection... "; flush stdout
#	writevsf "$name$ident.vsf" ; writevcf "$name$ident.vcf"
#	puts -nonewline "Output created, establishing link... "; flush stdout
#	for {set port 10000} { $port < 65000 } { incr port } {
#	    catch {imd connect $port} res
#	    if {$res == ""} break
#	}
#	if { $port==65000 } { puts "Failed." } else { puts "Done (now listening at port $port)." 
#	   # puts "    What you have to do now for a VMD connection:"
#	   # puts "    (1) Start vmd in current directory (best before running the script)."
#	   # puts "    (2) Enter on vmd command line: 'mol load psf $name$ident.psf pdb $name$ident.pdb'"
#	    set HOSTNAME [exec hostname]
#	   # puts "    (3) Enter on vmd command line: 'imd connect $HOSTNAME $port'"
#	   # puts "    (4) To have the chains coloured individually, set 'Coloring-Method' to 'ResName' in the 'Graphics'-menu"
#	    imd listen 0
#	    set vmdout_file [open "vmdoutput.script" "w"]
#	    puts $vmdout_file "mol load vsf $name$ident.vsf vcf $name$ident.vcf"
#	    puts $vmdout_file "rotate stop"
#	    puts $vmdout_file "imd connect $HOSTNAME $port"
#	    close $vmdout_file
#	    exec vmd -e vmdoutput.script &
#	}
#    }
#
#

#############################################################
#  Warmup Integration                                       #
#############################################################

puts "\nStart warmup integration:"
puts "At maximum $warm_n_times times $warm_steps steps"
#puts "Stop if minimal distance is larger than $min_dist"

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

#############################################################
#      Integration                                          #
#############################################################
#
#puts "Prepare integration: (tune parameters - wait...)"
#inter forcecap 0
#setmd time 0.0
## Set attractive LJ for monomer interactions (types 0 and 1
#for {set ia1 0} { $ia1 < 2 } { incr ia1 } {
#    for {set ia2 0} { $ia2 < 2 } { incr ia2 } {
#	inter $ia1 $ia2 lennard-jones $lj2_epsilon 1.0 $lj2_cut $lj2_shift 0.0
#    }
#}
##puts "[inter coulomb $bjerrum p3m tune accuracy $accuracy mesh 16]"
#
#puts "Interactions are now: {[inter]}"
#
##      Write blockfiles for restart
###########################################################
set trajectory [open ./data/m$l_polymers/$name$ident.config "w"]

blockfile $trajectory write variable {box_l time_step skin}
blockfile $trajectory write interactions
blockfile $trajectory write integrate
blockfile $trajectory write thermostat
flush $trajectory



# prepare observable output
set obs_file [open ./data/m$l_polymers/$name$ident.obs "w"]
puts $obs_file "\#$name$ident: Observables"
puts $obs_file "\#Time     R_E         R_G         E_TOT       E_KINT       E_POT"
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

    set f [open ./data/m$l_polymers/$name$ident\_dist.dat w]
    puts -nonewline $f "\#"
    puts $f "[analyze distribution 2 {0 1} 1.0 300.0 30 1 1]"
    flush $f
    close $f

##Coordinate
########################################################
    
    set vcf_file [open ./data/m$l_polymers/$name$ident\.vcf "w"]
    puts $vcf_file "\#Coordinate"
    writevcf $vcf_file 
    flush $vcf_file
    close $vcf_file

##Paraview
    
    set vtk_file ./data/m$l_polymers/$name$ident\.vtk
    #set vtk_file [open "$name$ident\.vtk" "w"]
    #puts $vtk_file "\#Partickle Visualization in paraview"
    writevtk $vtk_file 
    #flush $vtk_file
    #close $vtk_file

    set energy [analyze energy]
    set re_full [analyze re]
    set re [lindex $re_full 0]
    set rg [lindex [analyze rg] 0]

    puts $obs_file [format "%.3e %.5e %.5e %.5e %.5e %.5e" $time $re $rg [lindex [lindex $energy 0] 1] [lindex [lindex $energy 1] 1] [expr [analyze energy total] - [analyze energy kinetic]]]
    flush $obs_file

    if { $vmd_output=="yes" } { imd positions }
    #   write intermediate configuration
    if { $i%50==0 } {
	blockfile $trajectory write particles
	flush $trajectory
	incr j
    }
}

puts "\n Finished"
exit
