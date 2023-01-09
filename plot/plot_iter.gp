
set terminal png size 1200,700
set output "iter.png"

set xlabel "Nombre de points"

set ylabel "Nombre d'itérations"

set key left top

plot "../bin/ITER.dat" u 1:2 t "Richardson Alpha" with lines lc "blue", \
     "../bin/ITER.dat" u 1:3 t "Jacobi" with lines lc "web-green", \
     "../bin/ITER.dat" u 1:4 t "Gauß-Seidel" with lines lc "red"
