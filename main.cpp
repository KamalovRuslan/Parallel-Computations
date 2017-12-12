#include <exception>

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <string>

#include <mpi.h>
using namespace std;

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Poisson.h"

int main(int argc, char** argv){

	if (argc != 4) {
        cout << "Wrong arguments. Example: ./prac dimension eps out.txt" << endl;
        return 1;
    }

	int dimension;
	stringstream(argv[1]) >> dimension;

	double eps;
	stringstream(argv[2]) >> eps;

	string fout = string(argv[3]);

	try{
		Solver solver(dimension, eps, fout);
		solver.Solve(argc, argv);
	}
	catch (exception& e) {
        cout << e.what() << endl;
		return 1;
	}
	return 0;
}
