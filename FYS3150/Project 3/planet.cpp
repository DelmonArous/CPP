#include <iostream>
#include <cmath>
#include <string>
#include "planet.h"
#include "constants.h"

using namespace std;

// Constructor
planet::planet() {
}
planet::planet(string object, double dMass, double m_dist , double v_0, double angl) {
	name = object;
	mass = dMass;
	dist = m_dist/AU; // dim.less
	v0 = v_0/(AU/yr); // dim.less
	theta = angl;

	// Convert from polar to Cartesian coordinates
	xPos = cos(angl)*dist;
	yPos = sin(angl)*dist;
	Vx = -sin(angl)*v0;
	Vy = cos(angl)*v0;
}
