CXXFLAGS = -Wall -Werror -std=c++11 -g
fractals: fractals.o
	g++ $(CXXFLAGS) -o fractals fractals.o
fractals.o: fractals.cc
	g++ $(CXXFLAGS) -c fractals.cc
clean:
	rm  *.o fractals
