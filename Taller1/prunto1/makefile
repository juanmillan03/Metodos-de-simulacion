all: compile link gnu 

compile:
	g++ -O3 Modelosir_A.cpp -o Modelosir_A.x
	g++ Modelosir_C.cpp -o C.x

link: 
	Modelosir_A.x > datos.dat
	C.x

gnu:
	gnuplot  -p grafic.gp
	gnuplot  -p grafic3.gp
	rm -f *.o *.exe *.dat

clean:
	del/ -f *.x *.dat *.pdf