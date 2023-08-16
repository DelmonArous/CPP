#pragma once
#include "armadillo"
#include <string>

using namespace arma;
using namespace std;

class solver2D {
private:
	int xSteps, ySteps, tSteps; 
	double delta_x, delta_y, delta_t;
	string filename;
public:
	solver2D();
	void initialise();
	void explicit_scheme();
	void implicit_scheme();
	void analytic();
	void output(int, mat, string);
};