CC = g++
CFLAGS = -fopenmp

all: run

LBD3Q7tests: LBD3Q7tests.cpp
	$(CC) $(CFLAGS) -o $@ $<

run: LBD3Q7tests
	OMP_NUM_THREADS=4 ./LBD3Q7tests

clean:
	rm -f LBD3Q7tests *.dat

