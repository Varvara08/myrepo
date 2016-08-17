mol load vtf text.vtf
logfile vmd.log
mol new {/home/varya/mygit/UP/espresso/text.vtf} type {vtf} first 0 last -1 step 1 waitfor 1
mol modstyle 0 0 CPK 1.000000 0.300000 12.000000 12.000000
pbc box_draw
#draw color red
set box_length 7.0
draw triangle [list 0 $box_length 0] {0 0 0} [list $box_length 0 0]
draw triangle [list 0 $box_length 0] [list $box_length $box_length 0] [list $box_length 0 0]
#draw triangle {0 $box_length 0} {$box_length $box_length 0} {$box_length 0 0}

draw color red
draw sphere {3 3 3} radius 0.7
