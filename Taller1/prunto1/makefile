all: compile link gnu 

compile:
	g++ -O3 Modelosir_A.cpp -o A.exe
	g++ Modelosir_C.cpp -o C.exe

link: 
	./A.exe > datos.dat
	./C.x

gnu:
	gnuplot  -p grafic.gp
	gnuplot  -p grafic3.gp
	rm -f *.o *.exe *.dat

clean:
	rm -f  *.exe *.dat *.pdf