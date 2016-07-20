mol load vtf espresso_logo_online.vtf
imd connect localhost 10000

display resize 600 600
render options POV3 vmd_povray +W600 +H600 +ua -I%s +FN &
axes location off
color Display Background 8
# cup blue
color Name O 0
# saucer red
color Name N 1 
# steam silver
color Name S 6

animate goto 44
scale to 0.12
translate to 0 0.4 0.5

mol addrep 0

# cup and saucer
mol modselect 0 0 "not name S"
mol modstyle 0 0 CPK 3 0.3 8 6
mol modmaterial 0 0 Glossy

# steam
mol modselect 1 0 "name S"
mol modstyle 1 0 CPK 3 0.3 8 6
mol modmaterial 1 0 Glass2

