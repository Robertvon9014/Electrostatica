set title "Tiempo de ejecución vs Número de hilos"
set xlabel "Número de hilos"
set ylabel "Tiempo (s)"
set key top right
set grid
set xtics 1
set style data linespoints

plot \
    "SpeedUpM500.dat" using 1:2 title "M = 100" lt rgb "red", \
    "SpeedUpM1000.dat" using 1:2 title "M = 200" lt rgb "blue", \
    "SpeedUpM1500.dat" using 1:2 title "M = 300" lt rgb "green" \
    "SpeedUpM500.dat" using 1:(124.3248/$1) title "Ideal" lt rgb "black" dt 2

