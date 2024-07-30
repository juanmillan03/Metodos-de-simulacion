# Define una variable ARG por defecto vacía
Tmax ?=

all:
	g++ -O3 CA.cpp -o CA
	./CA $(Tmax)

Tasa:
	g++ -O3 CA.cpp -o CA
	./CA $(Tmax) > Resultados/Cuadricula_$(Tmax).dat 2> Resultados/sigma2.dat
Parallelo:
	g++ -O3 CA.cpp -o CA
	parallel './CA {} > Resultados/Cuadricula_{}.dat 2>/dev/null'  ::: 0 50 100 150 200 250 300 350 400 600 1000


