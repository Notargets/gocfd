# First plot to PNG
set terminal png size 1920,1080 font "Arial,14"
set output "plot_output.png"

# Set labels with specific fonts
set xlabel "X" font "Arial,16"
set ylabel "Value" font "Arial,16"
set title "Y, Mach, and Cp vs X" font "Arial,18"

# Set tics font
set xtics font "Arial,12"
set ytics font "Arial,12"

# Set key (legend) font
set key font "Arial,12"

set grid

# IMPORTANT: Set comma as the data separator
set datafile separator ","

# Define line styles with explicit dash patterns
set style line 1 lt 1 lw 2 lc rgb "red"      # Y (solid)
set style line 2 lt 1 lw 2 lc rgb "blue"     # Mach (solid)
set style line 3 lt 1 lw 2 lc rgb "green"    # Cp (solid)
#set style line 4 lt 1 dt (5,5) lw 2 lc rgb "blue"   # SU2 Mach (dashed: 5 on, 5 off)
#set style line 5 lt 1 dt (5,5) lw 2 lc rgb "green"  # SU2 Cp (dashed: 5 on, 5 off)
set style line 4 lt 0 lw 2 lc rgb "blue"     # SU2 Mach (dashed - uses built-in dash)
set style line 5 lt 0 lw 2 lc rgb "green"    # SU2 Cp (dashed - uses built-in dash)

# Plot both datasets with custom styles
plot 'plot.dat' using 1:2 with lines linestyle 1 title 'Y', \
     ''          using 1:3 with lines linestyle 2 title 'Mach', \
     ''          using 1:4 with lines linestyle 3 title 'Cp', \
     'RUN_SU2/mach_cp_data.dat' using 1:2 with lines linestyle 4 title 'SU2 Mach', \
     'RUN_SU2/mach_cp_data.dat' using 1:3 with lines linestyle 5 title 'SU2 Cp'

set terminal qt size 1920,1080 enhanced font 'Helvetica,10'
unset output

replot

pause -1
