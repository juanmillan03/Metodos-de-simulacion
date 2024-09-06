CC = g++
CFLAGS = -fopenmp

all: run

LBD3Q7tests: LBD3Q7tests.cpp
	$(CC) $(CFLAGS) -o $@ $<

run: LBD3Q7tests
	OMP_NUM_THREADS=4 ./LBD3Q7tests

clean:
	rm -f LBD3Q7tests *.dat
	rm -f D2/*.txt
	rm -f D3/*/*.txt
D2Q5:
	g++ -O3 LBOndas.cpp
	./a.out 
	python3 D2/gif.py
	eog D2/evolucion_sistema.gif
LBD3Q7_M: 
	g++ -O3 -fopenmp  LBD3Q7_millan.cpp -o LBM
	OMP_NUM_THREADS=16 ./LBM
	python3 D3/gif.py
