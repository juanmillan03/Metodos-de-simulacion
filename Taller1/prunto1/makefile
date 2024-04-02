all: py 

compile.x:
	g++ -O3 Modelosir_A.cpp -o Modelosir_A.x
	g++ -O3 Modelosir_C.cpp -o C.x

dat.x: compile.x
	Modelosir_A.x > datos.dat
	C.x

py:dat.x
	python grafica_C.py datos_C.dat
	python Modelosir_B.py Coronavirus_Colombia_y_Bogota.xlsx datos.dat
	del -f *.o *.x
