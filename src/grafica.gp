set title "Tiempo de ejecución vs Número de hilos"
set xlabel "Número de hilos"
set ylabel "Tiempo (s)"
set key top right
set grid
set xtics 1
set style data linespoints

plot \
    "tiempos_M500.dat" using 1:3 title "M = 100" lt rgb "red", \
    "tiempos_M1000.dat" using 1:3 title "M = 200" lt rgb "blue", \
    "tiempos_M1500.dat" using 1:3 title "M = 300" lt rgb "green"

