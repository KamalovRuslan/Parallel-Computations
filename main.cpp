#include <exception>

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include <mpi.h>
using namespace std;

#include "Poisson.h"

int main(int argc, char** argv){
	solver = Solver();
	try{
		solver.Solve();
	}
	catch (exception& e) {
        cout << e.what() << endl;
		return 1;
	}
	return 0;
}
