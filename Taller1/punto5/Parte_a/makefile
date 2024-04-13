inputs := 160 165 175 180 189 190 191
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
