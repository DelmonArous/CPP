#include <iostream>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <string>
#include "solver2D.h"
#include "constants.h"

solver2D::solver2D() {
}

void solver2D::initialise() {
	cout << "Delta x = "; cin >> delta_x;
	cout << "\nDelta t = "; cin >> delta_t;
	cout << "\nNumber of steps in t = "; cin >> tSteps;
	delta_y = delta_x;
	// Number of meshpoints in x and in y
	xSteps = (int) 1/delta_x;
	ySteps = xSteps;
}

void solver2D::explicit_scheme() {
	double x, y, t;
	double alpha = delta_t/(delta_x*delta_x);
	mat uOld(xSteps+1,ySteps+1), uNew(xSteps+1,ySteps+1);
	
	// Initial condition
	uOld.zeros(); uNew.zeros();
	for (int i = 0; i <= xSteps; i++) {
		x = i*delta_x;
		for (int j = 0; j <= ySteps; j++) {
			y = j*delta_y;
			uOld(i,j) = (1.-y)*exp(x);
		}
	}

	for (int k = 1; k <= tSteps; k++) { // time iteration
		t = k*delta_t;
		// Boundary conditions
		for (int i = 0; i <= xSteps; i++) {
			x = i*delta_x;
			uOld(i,0) = exp(x+t);
			uOld(i,ySteps) = 0.;
		}
		for (int j = 0; j <= ySteps; j++) {
			y = j*delta_y;
			uOld(0,j) = (1.-y)*exp(t);
			uOld(xSteps,j) = (1.-y)*exp(1.+t);
		}
		for (int i = 1; i < xSteps; i++) { // iterate x
			for (int j = 1; j < ySteps; j++) { // iterate y
				// Discretized diff eq.
				uNew(i,j) = alpha*(uOld(i-1,j) + uOld(i+1,j)+
					uOld(i,j-1) + uOld(i,j+1)) 
					+ (1-4*alpha)*uOld(i,j) ;
			}
		} // note that the boundaries are not changed
		// Replace the previous time solution with new		
		for (int i = 1; i < xSteps; i++) {
			for (int j = 1; j < ySteps; j++) {
				uOld(i,j) = uNew(i,j);
			}
		}
		// Write out to file
		if (k % 10 == 0) {
			filename = "Explicit_Iteration";
			filename += to_string(k) + ".txt";
			output(k, uOld, filename);
		}
	}
}

void solver2D::implicit_scheme() {
	double x, y, t, tOld;
	double alpha = delta_t/(delta_x*delta_x), max_error, acc_error, temp_error;
	int iteration, max_iter;
	mat uOld(xSteps+1,ySteps+1), uNew(xSteps+1,ySteps+1);
	
	max_iter = 10000; //Max number of allowed iteration
	acc_error = 0.00001; // Accepted error
/*	cout<<"\nMax number of iteration for Implicit method: "; cin>>max_iter;
	cout<<"\nMaximal Error treshold: "; cin>>acc_error;*/	

	// Initial condition
	uOld.zeros(); uNew.zeros();
	for (int i = 0; i <= xSteps; i++) {
		x = i*delta_x;
		for (int j = 0; j <= ySteps; j++) {
			y = j*delta_y;
			uOld(i,j) = (1.-y)*exp(x);
		}
	}

	for (int k = 1; k <= tSteps; k++) { // time iteration
		t = k*delta_t;
		
		// Boundary conditions
		for (int i = 0; i <= xSteps; i++) {
			x = i*delta_x;
			uOld(i,0) = exp(x+t);
			uOld(i,ySteps) = 0.;
		}
		for (int j = 0; j <= ySteps; j++) {
			y = j*delta_y;
			uOld(0,j) = (1.-y)*exp(t);
			uOld(xSteps,j) = (1.-y)*exp(1.+t);
		}

		uNew = uOld;// This is initial guess for timestep k
		iteration = 0;
		max_error = 1; // just to enter the while loop
		
		while( iteration < max_iter && max_error > acc_error){ // Jacobis method
			max_error = 0.0;
			temp_error = 0.0;
			for (int i = 1; i < xSteps; i++) { // iterate x
				for (int j = 1; j < ySteps; j++) { // iterate y
					// Discretized diff eq.
					uNew(i,j) = (1./(1+4*alpha))*(alpha*(uNew(i,j+1)+uNew(i,j-1)
						+uNew(i+1,j)+uNew(i-1,j))+uOld(i,j));
					temp_error = fabs(uNew(i,j)-uOld(i,j));
					if(max_error <= temp_error){max_error = temp_error;}

				}
			} 
			iteration++; 
			max_error /= pow(xSteps,2.); // Error acceptance esimate
		}

		// Write out to file
		if (k % 10 == 0) {
			cout<< iteration <<endl;
			filename = "Implicit_Iteration";
			filename += to_string(k) + ".txt";
			output(k, uNew, filename);
		}

		uOld = uNew;
	}
}

void solver2D::analytic() {
	double x, y, t;
	mat u(xSteps+1, ySteps+1);
	
	u.zeros();
	for (int k = 0; k <= tSteps; k++) { // time iteration
		t = k*delta_t;
		for (int j = 0; j <= ySteps; j++) { // iterate x
			y = j*delta_y;
			for (int i = 0; i <= xSteps; i++) { // iterate y
				x = i*delta_x;
				u(i,j) = (1. - y)*exp(x+t); // analytic solution
			}
		}
		if (k % 10 == 0) {
			filename = "analytic_t";
			filename += to_string(k) + ".txt";
			output(k, u, filename);
		}
	}
}

void solver2D::output(int step, mat matrix_2D, string filename) {
	ofstream outfile(filename);
	outfile << setiosflags(ios::showpoint | ios::uppercase);
	for (int j = 0; j <= ySteps; j++) {
		for (int i = 0; i <= xSteps; i++) {
			outfile << setw(15) << setprecision(8) << matrix_2D(i,j);
		}
		outfile << endl;
	}	
	outfile.close();
}