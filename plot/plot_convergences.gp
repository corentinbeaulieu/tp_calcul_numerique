
set terminal png size 1200,700
set output "convergences.png"

set xrange [0:100]
set xlabel "Nombre d'itérations"

set yrange [0:0.5]
set ylabel "Résidu relatif"


plot "bin/RESVEC_RICH.dat" t "Richardson" with lines, \
     "bin/RESVEC_JACOBI.dat" t "Jacobi" with lines lc "blue", \
     "bin/RESVEC_GAUSS_SEIDEL.dat" t "Gauß-Seidel" with lines lc "red"
