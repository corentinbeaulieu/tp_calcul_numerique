
set terminal png size 1200,700
set output "time_direct.png"

set xlabel "Nombre de points"

set ylabel "Temps (s)"

set logscale xy 2

set key left top

plot "bin/TIME_DIRECT.dat" u 1:2 t "dgbtrf+dgbtrs" with lines lc "blue", \
     "bin/TIME_DIRECT.dat" u 1:3 t "dgbsv" with lines lc "green", \
     "bin/TIME_DIRECT.dat" u 1:4 t "dgbtrftridiag" with lines lc "red"
