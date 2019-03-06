# Basic make file.

all: matrices

matrices: main.cpp matrices.cpp matrices.hpp
	g++ -g -Wall -Wextra -std=c++14 -o matrices main.cpp matrices.cpp