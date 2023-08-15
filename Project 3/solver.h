#pragma once

#include "planet.h"
#include "armadillo"

using namespace std;
using namespace arma;

class solver {
private:
	int n;
	double dt;
	vec uin;
public:
	planet *list;
	solver(int, double, vec);
	void addObject(int, planet);
	void initialize();
	vec forces(vec);
	vec RK4();
};

