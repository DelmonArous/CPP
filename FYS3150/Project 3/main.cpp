#include <stdio.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "armadillo"
#include "constants.h"
#include "planet.h"
#include "solver.h"

using namespace std;
using namespace arma;

int main() {

	int t_end, n, numberOfobjects;
	double t = 0.0, delta_t;
	cout << "Solving the differential equation of the n-body problem." << endl;
	cout << "Final time step in yr: "; cin >> t_end;
	cout << "Number of time/integration steps per year: "; cin >> n;
	cout << "Number of planets to include the system: "; cin >> numberOfobjects;
	numberOfobjects += 1; // include the Sun
	n *= t_end;
	delta_t = t_end/((double) n); // dimensionless timestep

	vec uin = zeros(4*numberOfobjects);
	planet *list = new planet[numberOfobjects];

	// How to call: planet(Name, Mass [kg], dist [m], v0 [m/s], angle [degrees])
	planet Sun("Sun", 1.99e30, 0.0, 0.0, 0.0);
	planet Mercury("Mercury", 3.30e23, 69.82e9, 38.86e3, 350.);
	planet Venus("Venus", 4.87e24, 108.94e9, 34.79e3, 45.);
	planet Earth("Earth", 5.97e24, AU, 29.29e3, 20.);
	planet Mars("Mars", 6.42e23, 249.23e9, 21.97e3, 250.);
	planet Jupiter("Jupiter", 1000*1.90e27, 816.62e9, 12.44e3, 285.);
	planet Saturn("Saturn", 5.68e26, 1514.50e9, 9.09e3, 0.);
	planet Uranus("Uranus", 8.68e25, 3003.62e9, 6.49e3, 75.);
	planet Neptune("Neptune", 1.02e26, 4545.67e9, 5.37e3, 295.);
	planet Pluto("Pluto", 1.31e22, 7375.93e9, 3.71e3, 295.);
	
	list[0] = Sun;
	list[1] = Mercury;
	list[2] = Venus;
	list[3] = Earth;
	list[4] = Mars;
	list[5] = Jupiter;
	list[6] = Saturn;
	list[7] = Uranus;
	list[8] = Neptune;
	list[9] = Pluto;
	
	solver test(numberOfobjects, delta_t, uin);

	for(int i = 0; i < numberOfobjects; i++) {
		test.addObject(i, list[i]);
	}
	test.initialize();

	fstream outfile;
    outfile.open("data_3body_m1000.txt", ios::out);

    while (t < t_end) {
        
		uin = test.RK4();

        t += delta_t;

        // Writes the position of all objects to file.
        for (int i = 0; i < numberOfobjects; i++) {
			outfile << uin(4*i+0) << "\t\t" << uin(4*i+1) << "\t\t";
        }
        outfile << endl;
    }
    outfile.close();

	return 0;
}