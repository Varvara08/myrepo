set term gif \
animate \
optimize \
delay 10 \
size 600, 600 \
background "#ffeedf" \
crop \
font "Times-Roman,14"
set output "animate.gif"
set size square
set xrange[0:200]
set yrange[0:0.2]
set xzeroaxis
set title "Distribution LB2_M=232"
do for [i=1:650] {
plot \
"CSBrush_m232_".i.".dat" u 1:7 lw 2 title sprintf("s = %d",i)
}
