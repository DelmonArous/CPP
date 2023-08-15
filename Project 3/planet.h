#pragma once
#include <string>

using namespace std;

class planet {
public:
	double xPos, yPos, dist, theta, v0, Vx, Vy, mass;
	string name;
	planet();
	planet(string, double, double, double, double);
};