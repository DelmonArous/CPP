#include <iostream>
#include <string>
#include "armadillo"
#include "MCsolver.h"
#include "solver2D.h"
#include "constants.h"
#include "lib.h"

using namespace arma;
using namespace std;


int main() {

	// The 1+1 dimension solver
	MCsolver test;
	test.initialise();
	test.MC_randomWalk();
	test.analytic_solution();

	// The 2+1 dimension solver
	solver2D solution2d = solver2D();
	solution2d.initialise();
	solution2d.analytic();
	solution2d.explicit_scheme();
	solution2d.implicit_scheme();



	return 0;
}