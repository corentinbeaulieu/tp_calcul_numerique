
set terminal png size 1200,700
set output "alpha.png"

set xlabel "alpha"

set ylabel "Number of iterations"

plot "bin/ALPHA.dat" notitle with lines
