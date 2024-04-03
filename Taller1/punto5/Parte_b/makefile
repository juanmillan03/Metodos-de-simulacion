inputs :=  157 165 190

all: compile link

compile:
	g++ -O3 Arena_B.cpp -o a.x

link:
	for n in $(inputs); do \
		./a.x "$$n" | gnuplot; \
	done
	rm -f *.x

clean:
	rm -f *.x