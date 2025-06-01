# First plot to PNG
set terminal png size 1920,1080 font "Arial,14"
set output "plot_output.png"


set title "Sod Shock Tube - Aligned Data"
set xlabel "X"
set ylabel "Solution Variables"
set key top right
set grid

plot '/tmp/tmp1' using 1:2 with lines lw 2 title 'Rho', \
     '' using 1:3 with lines lw 2 title 'RhoU', \
     '' using 1:4 with lines lw 2 title 'E', \
     '/tmp/tmp2' using 1:2 with points pt 7 ps 1.5 lc rgb 'red' title 'Analytic Rho', \
     '' using 1:3 with points pt 7 ps 1.5 lc rgb 'blue' title 'Analytic RhoU', \
     '' using 1:4 with points pt 7 ps 1.5 lc rgb 'green' title 'Analytic E'

set terminal qt size 1920,1080 enhanced font 'Helvetica,10'
unset output

replot

pause -1
