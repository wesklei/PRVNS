#!/bin/sh

CC=gcc

OBJECTS = mersenne.o 
SRCDIR = ./src-clstr

algorithm: $(OBJECTS)
	$(CC) $(OBJECTS) algorithm.c -Wall -Wextra -Wformat -Wformat-security -o algorithm -O3 -lm

mersenne.o: mersenne.c
	$(CC) mersenne.c -c -O3 -lm -Wall -Wextra

clean:
	rm algorithm  mersenne.o

