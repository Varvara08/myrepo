set term pngcairo dashed size 1280,960 font 'Times-Roman,18'
set output '| display png:-'

#set term eps size 6,4 font 'Times-Roman,18' 
#set output 'many_f_z2.eps'

#set title '{/Times-Italic f}=3'
set key right top

set yrange [0:]
set xrange [0:205]
#set xrange [0:30000]

set ylabel '{/Times-Italic V(z)=W(0,H)-W(z,H)}'
set ylabel '{/Times-Italic Terminal group}'
#set xlabel '{/Times-Italic z^2}'
set xlabel '{/Times-Italic z}'

set style line 2 lw 3.0 #lc rgb "#0B4C5F" #dt 2 pt 6 #
set style line 3 lw 1.5 #lc rgb "#F7D358" #dt 2 pt 6 #
set style line 5 lw 1.0 #lc rgb "#088A29" #dt 2 pt 6 #


p \
'sfbox_f3_100.dat' u 1:4 w l ls 2 dt 2 ti '{/Times-Italic f} = 3 {/Times-Italic {/Symbol s}} = 0.1',\
'sfbox_f3_200.dat' u 1:4 w l ls 2 dt 2 ti '0.2',\
'sfbox_f3_300.dat' u 1:4 w l ls 2 dt 2 ti '0.3',\
'check_sig100_z_pbp_pe_V1_np100.dat' u 1:3 w l ls 1 ti 'script  {/Times-Italic {/Symbol s}} = 0.1',\
'check_sig200_z_pbp_pe_V1_np100.dat' u 1:3 w l ls 1 ti 'script 0.2',\
'check_sig300_z_pbp_pe_V1_np100.dat' u 1:3 w l ls 1 ti 'script 0.3'
pause -1
#'sfbox_f4_100_man.dat' u ($1**2):2 w l ls 3 ti '{/Times-Italic f} = 4 0.1',\
#'sfbox_f4_200_man.dat' u ($1**2):2 w l ls 3 dt 3 ti '0.2',\
#'sfbox_f4_300_man.dat' u ($1**2):2 w l ls 3 dt 2 ti '0.3',\
#'sfbox_f8_100_man.dat' u ($1**2):2 w l ls 5 ti '{/Times-Italic f} = 8 0.1',\
#'sfbox_f8_200_man.dat' u ($1**2):2 w l ls 5 dt 3 ti '0.2'
#pause -1
