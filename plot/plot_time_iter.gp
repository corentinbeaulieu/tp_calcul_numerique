
set terminal png size 1200,700
set output "time_iter.png"

set xlabel "Nombre de points"
set xrange [14:]

set ylabel "Temps (s)"

set logscale xy 2

set key left top

plot "../bin/TIME_ITER.dat" u 1:2 t "Richardson Alpha" with lines lc "blue", \
     "../bin/TIME_ITER.dat" u 1:3 t "Jacobi" with lines lc "web-green", \
     "../bin/TIME_ITER.dat" u 1:4 t "Gauß-Seidel" with lines lc "red"
