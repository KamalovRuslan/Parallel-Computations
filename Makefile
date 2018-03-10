INCLUDES=$(wildcard *.h)
SOURCES=$(wildcard *.cpp)

.PHONY: clean images

compile:
	mpicxx main.cpp Poisson.cpp Poisson.h -o prac -Wall -O3

run:
	./prac
	python3 plotter.py
