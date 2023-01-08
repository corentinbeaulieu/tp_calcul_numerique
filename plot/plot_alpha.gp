
set terminal png size 1200,700
set output "alpha.png"

set xrange [0:0.52]
set xlabel "alpha"

set yrange [0:10100]
set ylabel "Nombre d'itérations"

plot "../bin/ALPHA.dat" notitle with lines
