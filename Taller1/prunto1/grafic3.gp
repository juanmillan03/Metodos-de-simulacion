reset session
set term pdf
set out "modelosir_C.pdf"
set xlabel "Tiempo(Días)"
set ylabel "Fraccion de Infectados activos "
set grid
set autoscale
set key right
set key nobox
#-------------------linea-punto-typeline-pointtype-
plot 'datos_C.dat' u 1:3 w l lt 7 title 'Simulación','Datos_reales.txt'u 1:2 w p lt 22 pt 22 ps 0.25  title 'Reales'