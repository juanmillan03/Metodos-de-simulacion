inputs := 160 180 190 195 200 205 210

all: compile link

compile:
	g++ -O3 Arena.cpp -o a.x

link:
	for n in $(inputs); do \
		./a.x "$$n" | gnuplot; \
	done
	rm -f *.x

clean:
	rm -f *.x
