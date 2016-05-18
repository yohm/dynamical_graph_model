CPP=g++
OPT=-O3 -Wall -DNDEBUG -std=c++11
# OPT=-O0 -Wall -pg
INCLUDE=

all: main.out

main.out: main.cpp dynamical_graph_model.cpp Makefile
	$(CPP) $(OPT) $(INCLUDE) main.cpp $(LIBS) -o $@

clean:
	rm -f *.out *~ *.bak
