#pragma once
#include "armadillo"
#include <string>

using namespace arma;
using namespace std;

class MCsolver {
private:
	int xSteps, tSteps, numberOfwalkers; 
	double delta_x, delta_t, move_probability, r;
	string filename;
public:
	MCsolver();
	void initialise();
	void analytic_solution();
	void MC_randomWalk();
	double gaussian_deviate(long *);
	void output_MC(vec, vec, vec, vec, int, long);
	void output(vec, string);
};

