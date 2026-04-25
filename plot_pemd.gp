# Usage:
# gnuplot -c plot_pemd.gp datafile outputfile logscale

datafile = ARG1
outfile  = ARG2
logscale = int(ARG3)

set terminal pngcairo size 1000,700
set output outfile

set xlabel 'k (a.u.)' font ', 15' enhanced
set ylabel 'P(k)' font ', 15' enhanced
set title 'Photoelectron Momentum Distribution (PEMD)' font ",16" enhanced

set format y "10^{%L}";

set grid
set size ratio 0.3182 1,1
set key left top
set xrange [-1.33:1.33]
set yrange [1.e-9 : 1.e-3]
set ytics  1.e-9, 1.e-2, 1.e-2
set xtics -1.0, .5, 1.0


if (logscale == 1) {
    set log y
#   plot datafile using 1:(log10($2)) with lines lw 2 lc 8 title 'P(k)'
    plot datafile using 1:($2) with lines lw 2 lc 8 title 'P(k)'
} else {
    plot datafile using 1:2 with lines lw 2 title 'P(k)'
}
