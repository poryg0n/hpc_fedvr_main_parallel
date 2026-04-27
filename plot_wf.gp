# Usage:
# gnuplot -c plot_pemd.gp datafile outputfile logscale

# --- ignore header lines
set datafile commentschars "#"

# --- allow external file input
#if (!exists("file")) file = "wf.dat"

datafile = ARG1
outfile  = ARG2
logscale = int(ARG3)

set terminal pngcairo size 1000,700

set output outfile

set grid
set size ratio 0.3182 1,1


#set title sprintf("Wavefunction: %s", file)
set xlabel "x"
set ylabel "Amplitude"
set grid
set key left top

set xrange [-5.33:5.33]


# columns:
# 1:x  2:Re  3:Im  4:|psi|^2  5:arg

plot \
    datafile using 1:2 w l lw 2 title "Re(ψ)", \
    datafile using 1:3 w l lw 2 title "Im(ψ)", \
    datafile using 1:4 w l lw 2 title "|ψ|^2"
