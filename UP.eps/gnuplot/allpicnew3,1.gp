set term gif \
	animate \
	optimize \
	delay 10 \
	size 1600,3200\
	background "#ffffff"\
	crop \
	font "Times-Roman,18"
set output "animate.gif"
#set autoscale
#set xrange[0:150]
#set yrange[0:0.2]
set xzeroaxis
#OPTIONS
########
f=2
q=f-1
n=50
h2n=25
########
#LOADING COLOR
#load '~/mygit/gnuplot-colorbrewer/sequential/Greens.plt'
#
###############
#set title "f=".f." n=".n." mn=2.06 var sig"
#set colorsequence default 
set style line 2 lw 4.0 #lc rgb "#0B4C5F" dt 2 pt 6 #
set style line 3 lw 4.0 #lc rgb "#F7D358" dt 2 pt 6 #
set style line 5 lw 4.0 #lc rgb "#088A29" dt 2 pt 6 #
#mypath="~/imc/StarBrush/Neutral/n".n."/q".q."den/chi0/mn2.06/"
mypath='~/imc/StarBrush/Neutral/Poty/f'.f.'/h2n0.'.h2n.'/'
#mypath1='~/imc/StarBrush/Neutral/Poty/f'.f.'/h2n0.'.h2n.'/from_sfbox/'
#mypath1='~/imc/StarBrush/Neutral/Poty/f'.f.'/h2n0.'.h2n.'/M254/'
#mypath2='~/imc/StarBrush/Neutral/Poty/f'.f.'/h2n0.'.h2n.'/M256/'
mypath='~/imc/StarBrush/Neutral/Poty/f'.f.'/h2n0.'.h2n.'/M435/'
do for [i=100:550] {
set multiplot layout 3, 1 title 'Line chain in brush f=2_{254,255,256}'
set tmargin 2
#set title "{/Symbol \152}(z)"
set key
set xrange [0:80]
set yrange [0:18]
set ylabel 'Density profile linear chain'
set xlabel 'z'
plot mypath.'h2n0.'.h2n.'_f'.f.'8_z_phi_pe_m'.i.'.dat' u 1:2 w l ls 2 title sprintf ("M_{435} = %d",i)
#
set xrange [0:80]
set yrange [0:0.15]
set ylabel 'Chain`s end distribution'
set xlabel 'z'
#set title "Pe(z)"
set key
plot mypath.'h2n0.'.h2n.'_f'.f.'8_z_phi_pe_m'.i.'.dat' u 1:3 w l ls 2 title sprintf ("M_{435} = %d",i),\
#
#set title "V(z)"
set key
set xrange [0:80]
set yrange [0:0.3]
set ylabel 'Potention, V(0.'.h2n.')'
set xlabel 'z'
plot mypath.'h2n0.'.h2n.'_f'.f.'8_z_phi_pe_m'.i.'.dat' u 1:4 w l ls 2 title sprintf ("M_{435} = %d",i),\
}
