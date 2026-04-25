# Usage:
# gnuplot -c plot_hhg.gp datafile outputfile mode

# mode:
# 1 = dipole
# 2 = momentum
# 3 = acceleration
# 4 = IBP
# 5 = all

datafile = ARG1
outfile  = ARG2
mode     = int(ARG3)

set format y "10^{%L}";

set terminal pngcairo size 1000,700
set output outfile

set xlabel '{/Symbol w} (a.u.)' font ', 16'
set ylabel 'Q^{(c)}({/Symbol w})' font ', 16'
set title 'HHG Spectrum'

set logscale y
set grid
set key left bottom


set size ratio  0.3182 1,1
set xrange [.5:1.13]
#set yrange [1.e-14:1.e-2]
set xtics .5, .1, 1.1
set ytics 1.e-13, 1.e-2, 1.e-3

#if (mode == 1) {
#    plot datafile using 1:2 with lines lw 2 lc 8 title 'Eq.(104)'
#} if (mode == 2) {
#    plot datafile using 1:3 with lines lw 2 title 'Eq.(80)'
#} if (mode == 3) {
#    plot datafile using 1:4 with lines lw 2 title 'Eq.(40)'
#} if (mode == 4) {
#    plot datafile using 1:5 with lines lw 2 title 'IBP'
#} if (mode == 5) {
#    plot datafile using 1:2 with lines lw 2 title 'Eq.(104)', \
#         datafile using 1:3 with lines lw 2 title 'Eq.(80)', \
#         datafile using 1:4 with lines lw 2 title 'Eq.(1)'
#}

    plot datafile using 1:2 with lines lw 2 title 'Eq.(104)', \
         datafile using 1:3 with lines lw 2 title 'Eq.(80)', \
         datafile using 1:4 with lines lw 2 title 'Eq.(1)', \
#        datafile using 1:5 with lines lw 2 title 'IBP'
