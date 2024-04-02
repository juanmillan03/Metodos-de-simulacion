all: compile link

compile:
	g++ -O3 Arena.cpp -o a.x
link: 
	a.x | gnuplot
clean :
	rm -f *.x

# for name in Alice Bob Charlie
# do
#     echo "Hello, $name"
# done