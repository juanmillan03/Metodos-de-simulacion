reset session
set term pdf
set out "modelosir_A.pdf"
set xlabel "Tiempo(Días)"
set ylabel "Susceptibles, Infectados and Recuperados "
set grid
set autoscale
set title "SIR epidemic without births or deaths"
set key center
set key right
#-------------------linea-punto-typeline-pointtype-
plot 'datos.dat'u 1:2 w lp lt 14 pt 7 ps 0.5 title 'Susceptibles', \
     ''u 1:3 w lp lt 7 pt 7 ps 0.5 title 'Infectados', \
     ''u 1:4 w lp lt 2 pt 7 ps 0.5 title 'Recuperados' \