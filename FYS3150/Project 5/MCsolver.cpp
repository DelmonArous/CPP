#include <iostream>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <string>
#include "armadillo"
#include "MCsolver.h"
#include "constants.h"
#include "lib.h"

using namespace arma;
using namespace std;

MCsolver::MCsolver() {
}

void MCsolver::initialise() {
	cout << "Delta t = "; cin >> delta_t;
	//delta_t = delta_x*delta_x/2.;
	delta_x = sqrt(2.*delta_t); // D = 1 = delta_x^2/2*delta_t
	xSteps = (int) 1/delta_x; // Number of meshpoints in x
	cout << "\nNumber of steps in t = "; cin >> tSteps;
	cout << "\nNumber of random walkers = "; cin >> numberOfwalkers;
	cout << "\nMove probability = "; cin >> move_probability;
}

void MCsolver::analytic_solution() {
	double x, t;
	vec u(xSteps+1);
	u.zeros(); // initial condition
	// Boundary conditions
	u(0) = 1.; u(xSteps) = 0.;

	for (int k = 1; k <= tSteps ; k++) { // time iteration
		t = k*delta_t;
		for (int i = 1; i < xSteps; i++) { // iterate x
			x = i*delta_x;
			double sum = 0.;
			for (int n = 1; n < 100; n++) {
				sum += (2./(n*pi))*sin(n*pi*x)*pow(e, -n*n*pi*pi*t);
			}
			u(i) = 1. - x - sum;
		}
		// Print out to file
		if (k % 10 == 0) {
			filename = "analytic_t";
			filename += to_string(k) + "_deltat" + to_string(delta_t) + ".txt";
			output(u, filename); 
		}
	}
}

void MCsolver::MC_randomWalk() {
	long idum = -1; // initialise random number generator
	for (int seed_nr = 1; seed_nr <= 3; seed_nr++) {
		vec x(numberOfwalkers), walk_cumulative(tSteps+1), walk_cumulative2(tSteps+1);
		vec xAverage(tSteps+1), x2Average(tSteps+1), variance(tSteps+1);
		vec probHistogram(xSteps+1), uHistogram(xSteps+1);
		// Initialization of vectors
		x.zeros(); walk_cumulative.zeros(); walk_cumulative2.zeros();
		xAverage.zeros(); x2Average.zeros(); variance.zeros();
		probHistogram.zeros(); uHistogram.zeros();

		long seed_value = idum;
		int counter_x0, counter_x1, counter_zeros;
		for (int step = 1; step <= tSteps; step++) { // loop over all time steps 
			counter_x0 = 0; counter_x1 = 0; counter_zeros = 0;
			for (int walker = 0; walker < (int) x.n_elem; walker++) { // loop over number of walkers
				double r = ran2(&idum); // obtain a floating number r between [0,1]
				double l = delta_x*gaussian_deviate(&idum); // The jump-step
				if (r <= move_probability) {x(walker) += l;}
				else {x(walker) -= l;}

				if (x(walker) < 0.) {counter_x0++;}
				// Problem when the particles are in position x=l and then moves -l
				else if (x(walker) == 0. || x(walker) < pow(10,-15)) {counter_zeros++;}
				else if (x(walker) >= 1.) {counter_x1++;}
			}
			// This is our model of a random walker
			x = sort(x, 0);
			x.resize(x.n_elem - counter_x1); // remove particles with positions x >= 1.0
			x = sort(x, 1);
			x.resize(x.n_elem - counter_x0); // remove particles with positions x < 0.0
			x.resize(x.n_elem - counter_zeros); // remove particles with position x = 0.0
			// Add particles to position x = 0.0 since the density is constant here 
			x.resize(x.n_elem + numberOfwalkers);

			for (int walker = 0; walker < (int) x.n_elem; walker++) {
			walk_cumulative(step) += x(walker);
			walk_cumulative2(step) += x(walker)*x(walker);
			}
			// Average the displacements and their variances
			xAverage(step) = walk_cumulative(step)/((double) x.n_elem);
			x2Average(step) = walk_cumulative2(step)/((double) x.n_elem);
			variance(step) = x2Average(step) - xAverage(step)*xAverage(step);
		}
		// Estimate the probability by counting number
		// number walkers at every position on the final time step
		double norm = 0.;
		for (int i = 0; i < (int) probHistogram.n_elem ; i++) {
			int counter = 0;
			for (int walker = 0;  walker < (int) x.n_elem; walker++) {
				if (i/((double) probHistogram.n_elem-1.) <= x(walker) && x(walker) < (i+1.)/((double) probHistogram.n_elem-1.)) {
					counter++;
				}
			}
			probHistogram(i) = counter;
			norm += probHistogram(i);
		} 
		uHistogram = probHistogram/probHistogram(0);
		probHistogram = probHistogram/norm;
		output_MC(xAverage, variance, probHistogram, uHistogram, seed_nr, seed_value); // write out to file
	}
}

double MCsolver::gaussian_deviate(long * idum) {
  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;

  if ( idum < 0) iset = 0;
  if (iset == 0) {
	  do {
      v1 = 2.*ran2(idum) -1.0;
      v2 = 2.*ran2(idum) -1.0;
      rsq = v1*v1+v2*v2;
	  } while (rsq >= 1.0 || rsq == 0.);
	  fac = sqrt(-2.*log(rsq)/rsq);
	  gset = v1*fac;
	  iset = 1;
	  return v2*fac;
  } 
  else {
	  iset =0;
	  return gset;
  }
} // end function for gaussian deviates

void MCsolver::output_MC(vec xAverage, vec variance, vec w, vec histogram, int seed, long idum) {
	filename = "xAverage_variance_tSteps" + to_string(tSteps) + 
		"_deltat" + to_string(delta_t) + "_seed" + to_string(seed) + ".txt";
	ofstream ofile(filename);
	for (int k = 1; k <= tSteps; k++) {
		ofile << setiosflags(ios::showpoint | ios::uppercase);
		ofile << setw(6) << k; 
		ofile << setw(15) << setprecision(8) << xAverage(k);
		ofile << setw(15) << setprecision(8) << variance(k) << endl;
	}
	ofile.close();

	filename = "probability_tSteps" + to_string(tSteps) + 
		"_deltat" + to_string(delta_t) + "_seed" + to_string(seed) + ".txt";
	ofstream probfile(filename);
	for (int i = 0; i < (int) w.n_elem; i++) {
		probfile << setiosflags(ios::showpoint | ios::uppercase);
		probfile << setw(6) << i; 
		probfile << setw(15) << setprecision(8) << w(i) << endl;
	}
	probfile.close();

	filename = "histogram_tSteps" + to_string(tSteps) + 
		"_deltat" + to_string(delta_t) + "_seed" + to_string(seed) + ".txt";
	ofstream histfile(filename);
	for (int i = 0; i < (int) histogram.n_elem; i++) {
		histfile << setiosflags(ios::showpoint | ios::uppercase);
		histfile << setw(6) << i; 
		histfile << setw(15) << setprecision(8) << histogram(i) << endl;
	}
	histfile.close();

	filename = "seed__gaussReal_tSteps" + to_string(tSteps) + "_deltat" +
		to_string(delta_t) + ".txt";
	ofstream seedfile(filename, ios::app);
	seedfile << setiosflags(ios::showpoint | ios::uppercase);
	seedfile << setw(15) << setprecision(8) << idum << endl;
	seedfile.close();
}

void MCsolver::output(vec u, string filename) {
	ofstream outfile(filename);
	outfile << setiosflags (ios::showpoint | ios::uppercase);
	for (int i = 0; i <= xSteps; i++) {
		outfile << setw(15) << setprecision(8) << i*delta_x;
		outfile << setw(15) << setprecision(8) << u(i) << endl;
	}
	outfile.close();
}