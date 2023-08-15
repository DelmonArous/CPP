#include <cmath>
#include <iostream>
#include "solver.h"
#include "planet.h"
#include "constants.h"
#include "armadillo"

using namespace std;
using namespace arma;

solver::solver(int numberOfobjects, double timeStep, vec u) {
	n = numberOfobjects;
	dt = timeStep;
	list = new planet[n];
	uin = u;
}

void solver::addObject(int m, planet object) {
	list[m] = object;
}

void solver::initialize() {
	int k = 4;
	// Initializing the vector u(x1, y1, v1_x, v1_y, x2, y2, v2_x, v2_y, ..., xn, yn, vn_x, vn_y)
	for (int i = 0; i < n; i++) {
        this->uin(i*k + 0) = this->list[i].xPos;
        this->uin(i*k + 1) = this->list[i].yPos;
        this->uin(i*k + 2) = this->list[i].Vx;
        this->uin(i*k + 3) = this->list[i].Vy;
    }
    // cout << uin << endl;
}

vec solver::forces(vec uin) {
	vec Forces = zeros(2*n);
    vec dudt = zeros(4*n);
	double m_sun = list[0].mass;
    // Compute r and F between bodies.
    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
				double x = uin[4*i] - uin[4*j];
				double y = uin[4*i+1] - uin[4*j+1];
				double r = sqrt(x*x + y*y);
				double f = -(G * list[i].mass * list[j].mass) / (m_sun * r * r * r);
            
				Forces(2*i+0) += f*x;		// x-component of the force
				Forces(2*i+1) += f*y;      // y-component of the force
				Forces(2*j+0) -= Forces(2*i+0);    // Newtons third law
				Forces(2*j+1) -= Forces(2*i+1);    // Newtons third law
		}
	}

    for (int i = 0; i < n; i++) {
        double m = list[i].mass;
		// Compute dv_x and dv_y
        dudt(4*i+2) = Forces(2*i+0) / m; // dv_x/dt
		dudt(4*i+3) = Forces(2*i+1) / m; // dv_y/dt
		// Compute dx and dy
        dudt(4*i+0) = uin(4*i+2); // dx/dt
        dudt(4*i+1) = uin(4*i+3); // dy/dt
	}
	return dudt;
}

vec solver::RK4() {
	vec k1(4*n), k2(4*n), k3(4*n), k4(4*n);
    k1 = forces(uin) * dt;
    k2 = forces(uin + 0.5 * k1) * dt;
    k3 = forces(uin + 0.5 * k2) * dt;
    k4 = forces(uin + k3) * dt;
    uin += (1.0/6) * (k1 + 2 * (k2 + k3) + k4);
	return uin;
}