CC=g++
CFLAGS=-O2 -ansi -pedantic -Wall

EXE=Ab_lattice

all: $(EXE)

Ab_lattice: Ab_lattice.cpp field.hpp vec2.hpp
	$(CC) -o Ab_lattice Ab_lattice.cpp $(CFLAGS)
	chmod +x Ab_lattice

clean:
	rm $(EXE) 2>&-
