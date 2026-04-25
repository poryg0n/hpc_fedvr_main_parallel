set terminal pngcairo size 1000,700

datafile = ARG1
outfile  = ARG2
#mode     = int(ARG3)

set terminal pngcairo size 1000,700
set output outfile

set output outfile
set title "p_{k0} consistency check"
set xlabel "k (a.u.)" font ", 15"
set ylabel "p_{k0}" font ", 15"

set size ratio  0.3182 1,1
#set xrange [.5:1.13]
#set yrange [1.e-14:1.e-2]
#set xtics .5, .1, 1.1
#set ytics 1.e-13, 1.e-2, 1.e-3


set grid
set key left top

plot \
datafile using 1:2 with lines lw 2 title "Re(p_{k0})", \
datafile using 1:3 with lines lw 2 title "Re(p_{0k})", \
datafile using 1:4 with lines lw 2 title "Re(d_{k0} -> p)", \
datafile using 1:5 with lines dt 2 lc 8 title "Eq.(102a)"
