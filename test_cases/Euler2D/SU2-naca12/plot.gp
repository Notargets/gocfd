# Set terminal with size and font
set terminal qt size 1920,1080 font "Arial,14"

# Set labels with specific fonts
set xlabel "X" font "Arial,16"
set ylabel "Value" font "Arial,16"
set title "Y, Mach, and Cp vs X" font "Arial,18"

# Set tics font
set xtics font "Arial,12"
set ytics font "Arial,12"

# Set key (legend) font
set key font "Arial,12"

set xlabel "X"
set ylabel "Value"
set title "Y, Mach, and Cp vs X"
set grid

plot 'plot.dat' using 1:2 with lines title 'Y', \
     ''             using 1:3 with lines title 'Mach', \
     ''             using 1:4 with lines title 'Cp'

pause -1
