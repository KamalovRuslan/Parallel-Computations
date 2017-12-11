#include <exception>

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include <mpi.h>
using namespace std;

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Poisson.h"

int main(int argc, char** argv){
	try{
		Solver solver;
		solver.Solve(argc, argv);
	}
	catch (exception& e) {
        cout << e.what() << endl;
		return 1;
	}
	return 0;
}
