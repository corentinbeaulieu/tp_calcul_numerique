
set terminal png size 1200,700
set output "convergences.png"

set xrange [0:125]
set xlabel "Number of iterations"

set yrange [0:0.5]
set ylabel "Forward Error"


plot "bin/RESVEC_RICH.dat" t "Richardson" with lines, \
     "bin/RESVEC_JACOBI.dat" t "Jacobi" with lines, \
     "bin/RESVEC_GAUSS_SEIDEL.dat" t "Gauß-Seidel" with lines
