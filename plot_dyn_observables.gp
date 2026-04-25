# Usage:
# gnuplot -c plot_density.gp datafile outputfile

datafile = ARG1
outfile  = ARG2
#mode     = int(ARG3)

#set format y "10^{%L}";

set terminal pngcairo size 1000,700
set output outfile

#set logscale y
#set grid
#set key left bottom

set size ratio  0.3182 1,1
set xrange [-1120:1120]
#set yrange [-3.e-1:3.e-1]
#set xtics .5, .1, 1.1
#set ytics 1.e-13, 1.e-2, 1.e-3

set grid

set key top right

set xlabel 'time (a.u.)'
set ylabel 'd(t)'
set grid

set key top right


#set title "Observables vs Time"
set xlabel "time"
set grid

set multiplot layout 3,1

# --- Energy ---
set ylabel "Energy"
plot datafile using 1:2 w l lw 2 title "Energy"

# --- Dipole ---
set ylabel "Dipole"
plot datafile using 1:3 w l lw 2 title "Re[d(t)]", \
     datafile using 1:5 w l title "Im[d(t)]"

# --- Momentum ---
set ylabel "Momentum"
plot datafile using 1:4 w l lw 2 title "Re[p(t)]", \
     datafile using 1:6 w l title "Im[p(t)]"

unset multiplot
