#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "armadillo"

using namespace arma;
using namespace std;

void initialize(int &xSteps, int &tSteps, double &delta_x, double &delta_t);
void explicit_scheme(int xSteps, int tSteps, double delta_x, double delta_t);
void implicit_scheme(int xSteps, int tSteps, double delta_x, double delta_t);
void Crank_Nicolson_scheme(int xSteps, int tSteps, double delta_x, double delta_t);
void tridiag(vec a, vec b, vec c, vec &Vold, vec &Vnew, int xSteps);
void analytic_solution(int xSteps, int tSteps, double delta_x, double delta_t);
void output(int counter, int xSteps, int tSteps, double delta_x, double delta_t, vec Vnew, int meth_nr);
double func(double x);

int main() {
	int xSteps, tSteps;
	double delta_x, delta_t;

	// Read in data
	initialize(xSteps, tSteps, delta_x, delta_t);
	// Run forward Euler scheme
	explicit_scheme(xSteps, tSteps, delta_x, delta_t);
	// Run backward Euler scheme
	implicit_scheme(xSteps, tSteps, delta_x, delta_t);
	// Run Crank-Nicolson scheme
	Crank_Nicolson_scheme(xSteps, tSteps, delta_x, delta_t);
	// Simulate the closed-form soluton
	analytic_solution(xSteps, tSteps, delta_x, delta_t);
	return 0;
}

void initialize(int &xSteps, int &tSteps, double &delta_x, double &delta_t) {
	cout << "Enter Delta x: "; cin >> delta_x;
	cout << "\nEnter Delta t: "; cin >> delta_t;
	cout << "\nEnter number of steps in t: "; cin >> tSteps;
	xSteps = 1/delta_x; // Number of mesh points in x
}

void explicit_scheme(int xSteps, int tSteps, double delta_x, double delta_t) {
	double x, alpha = delta_t/(delta_x*delta_x);
	int meth_nr = 1;
	vec Vold(xSteps+1), Vnew(xSteps+1), u(xSteps+1); // n+1 is the number of mesh points in x
	// Boundary conditions
	Vold(0) = 0.; Vold(xSteps) = 0.;
	Vnew(0) = 0.; Vnew(xSteps) = 0.;
	u.zeros();
	u(0) = 1.; u(xSteps) = 0.;

	// Initial condition
	for (int i = 1; i < xSteps; i++) {
		x = i*delta_x;
		Vold(i) = func(x);
	}
	int counter = 0;
	output(counter, xSteps, tSteps, delta_x, delta_t, u, meth_nr);

	// Time iteration
	for (int j = 1; j <= tSteps; j++) {
		for (int i = 1; i < xSteps; i++) {
			// Discretized diff eq
			Vnew(i) = alpha*Vold(i-1) + (1 - 2*alpha)*Vold(i) + alpha*Vold(i+1);
		} // note that the boundaries are not changed.
		// replace previous time solution with new
		for (int i = 0; i <= xSteps; i++) {
			x = i*delta_x;
			u(i) = 1.-x + Vnew(i);
			Vold(i) = Vnew(i);
		}
		// Print out to file
		output(j, xSteps, tSteps, delta_x, delta_t, u,meth_nr);
	}
}

void implicit_scheme(int xSteps, int tSteps, double delta_x, double delta_t) {
	double x, alpha = delta_t/(delta_x*delta_x);
	int meth_nr = 2;
	vec Vold(xSteps+1), Vnew(xSteps+1), u(xSteps+1);
	// Boundary condition
	Vold(0) = 0.; Vold(xSteps) = 0.; 
	Vnew(0) = 0.; Vnew(xSteps) = 0.;
	u.zeros();
	u(0) = 1.; u(xSteps) = 0.;

	// Initialize vectors
	for (int i = 1; i < xSteps; i++) {
		x = i*delta_x;
		Vold(i) = func(x); // initial condition
		Vnew(i) = Vold(i);
		//u(i) = 1.-x + Vnew(i);
	}
	int counter = 0;
	output(counter, xSteps, tSteps, delta_x, delta_t, u, meth_nr);

	// Setting up the matrix A, only constants
	vec a(xSteps), b(xSteps+1), c(xSteps);
	a.fill(-alpha); c.fill(-alpha); // a and c contains the off-diagonal elements in matrix A
	b.fill(1 + 2*alpha); // b constains the diagonal elements in matrix A

	// Time iteration
	for (int j = 1; j <= tSteps; j++) {
		// Here we solve the tridiagonal linear set of equations
		tridiag(a, b, c, Vold, Vnew, xSteps);
		Vnew(0) = 0.; Vnew(xSteps) = 0.; // boundary conditions
		// replace previous time solution with new
		for (int i = 0; i <= xSteps; i++) {
			x = i*delta_x;
			u(i) = 1.-x + Vnew(i);
			Vold(i) = Vnew(i);
		}
		// Print out to file
		output(j, xSteps, tSteps, delta_x, delta_t, u,meth_nr);
	}
}

void Crank_Nicolson_scheme(int xSteps, int tSteps, double delta_x, double delta_t) {
	double x, alpha = delta_t/(delta_x*delta_x);
	int meth_nr = 3;
	vec Vold(xSteps+1), Vnew(xSteps+1), u(xSteps+1);
	// Boundary condition
	Vnew(0) = 0.; Vnew(xSteps) = 0.;
	u.zeros();
	u(0) = 1.; u(xSteps) = 0.;

	// Initialize vector
	for (int i = 1; i < xSteps; i++) {
		x = i*delta_x;
		Vnew(i) = func(x); // initial condition
	}
	int counter = 0;
	output(counter, xSteps, tSteps, delta_x, delta_t, u, meth_nr);

	// Setting up the matrix A = 2*I + alpha*B
	vec a(xSteps), b(xSteps+1), c(xSteps);
	a.fill(-alpha); c.fill(-alpha);// a and c contains the off-diagonal elements in matrix A
	b.fill(2 + 2*alpha); // b contains the diagonal elements in matrix A

	// Time iteration
	for (int j = 1; j <= tSteps; j++) {
		// Calculate r for use in tridag, right hand side of the Crank Nicolson method
		for (int i = 1; i < xSteps; i++) {
			Vold(i) = alpha*Vnew(i-1) + (2 - 2*alpha)*Vnew(i) + alpha*Vnew(i+1);
		}
		Vold(0) = 0.; Vold(xSteps) = 0.; // boundary condition
		// Then solve the tridiagonal matrix
		tridiag(a, b, c, Vold, Vnew, xSteps);
		Vnew(0) = 0.; Vnew(xSteps) = 0.; // boundary condition
		// replace previous time solution with new
		for (int i = 0; i <= xSteps; i++) {
			x = i*delta_x;
			u(i) = 1.-x + Vnew(i);
			Vold(i) = Vnew(i);
		}
		// Print out to file
		output(j, xSteps, tSteps, delta_x, delta_t, u,meth_nr);
	}
}

void tridiag(vec a, vec b, vec c, vec &Vold, vec &Vnew, int xSteps) {
	double m;
	// Forward substitution
	for (int i = 1 ; i < xSteps; i++) {
		m = a(i)/b(i-1);
		b(i) -= m*c(i-1);
		Vold(i) -= m*Vold(i-1);
	}
	// Backward substitution
	Vnew(xSteps) = Vold(xSteps)/b(xSteps);
	for (int i = xSteps-1; i > 0 ; i--) {
		Vnew(i) = (Vold(i) - c(i)*Vnew(i+1))/b(i);
	}
}

void analytic_solution(int xSteps, int tSteps, double delta_x, double delta_t) {
	double x, t;
	int meth_nr = 0;
	const double pi = acos(-1.);
	const double e = 2.71828183;
	vec v(xSteps+1), u(xSteps+1);
	// Boundary conditions
	v(0) = 0.; v(xSteps) = 0.;
	u.zeros();
	u(0) = 1.; u(xSteps) = 0.;
	for (int i = 1; i < xSteps; i++) {
		x = i*delta_x;
		v(i) = x-1.; // initial condition
	}
	int counter = 0;
	output(counter, xSteps, tSteps, delta_x, delta_t, u,meth_nr);

	// Time iteration
	double sum;
	for (int j = 1; j <= tSteps; j++) {
		t = j*delta_t;
		for (int i = 1; i < xSteps; i++) {
			x = i*delta_x;
			sum = 0.;
			for (int n = 1; n < 100; n++) {
				sum += (2/(n*pi))*sin(n*pi*x)*pow(e, -n*n*pi*pi*t);
			}
			u(i) = 1. - x - sum;
		}
		// Print out to file
		output(j, xSteps, tSteps, delta_x, delta_t, u,meth_nr);
	}
}

void output(int counter, int xSteps, int tSteps, double delta_x, double delta_t, vec Vnew, int meth_nr) {
	string filename;
	ofstream outfile;
	if(meth_nr==0){
		if (counter == 0) {
			/*cout << "\nEnter the output filename: "; cin >> filename;
			filename += ".txt";*/
			outfile.open("analytic_test.txt", ios::out); // filename.c_str()
			outfile << setiosflags(ios::showpoint|ios::uppercase);
			for (int i = 0; i <= xSteps; i++) {
				outfile << setw(15) << setprecision(8) << i*delta_x << "\t\t";
			}
			outfile << endl;
		}
		else {
			outfile.open("analytic_test.txt", ios::app);
		}
		for (int i = 0; i <= xSteps; i++) {
			outfile << setw(15) << setprecision(8) << Vnew(i) << "\t\t";
		}
		outfile << endl;
		outfile.close();
	}
		if(meth_nr==1){
		if (counter == 0) {
			/*cout << "\nEnter the output filename: "; cin >> filename;
			filename += ".txt";*/
			outfile.open("explicit_test.txt", ios::out); // filename.c_str()
			outfile << setiosflags(ios::showpoint|ios::uppercase);
			for (int i = 0; i <= xSteps; i++) {
				outfile << setw(15) << setprecision(8) << i*delta_x << "\t\t";
			}
			outfile << endl;
		}
		else {
			outfile.open("explicit_test.txt", ios::app);
		}
		for (int i = 0; i <= xSteps; i++) {
			outfile << setw(15) << setprecision(8) << Vnew(i) << "\t\t";
		}
		outfile << endl;
		outfile.close();
	}
		if(meth_nr==2){
			if (counter == 0) {
				/*cout << "\nEnter the output filename: "; cin >> filename;
				filename += ".txt";*/
				outfile.open("implicit_test.txt", ios::out); // filename.c_str()
				outfile << setiosflags(ios::showpoint|ios::uppercase);
				for (int i = 0; i <= xSteps; i++) {
				outfile << setw(15) << setprecision(8) << i*delta_x << "\t\t";
				}
				outfile << endl;
			}
			else {
				outfile.open("implicit_test.txt", ios::app);
			}
			for (int i = 0; i <= xSteps; i++) {
				outfile << setw(15) << setprecision(8) << Vnew(i) << "\t\t";
			}
			outfile << endl;
			outfile.close();
		}
		else
			if (counter == 0) {
				/*cout << "\nEnter the output filename: "; cin >> filename;
				filename += ".txt";*/
				outfile.open("CN_test.txt", ios::out); // filename.c_str()
				outfile << setiosflags(ios::showpoint|ios::uppercase);
				for (int i = 0; i <= xSteps; i++) {
				outfile << setw(15) << setprecision(8) << i*delta_x << "\t\t";
				}
				outfile << endl;
			}
			else {
				outfile.open("CN_test.txt", ios::app);
			}
			for (int i = 0; i <= xSteps; i++) {
				outfile << setw(15) << setprecision(8) << Vnew(i) << "\t\t";
			}
			outfile << endl;
			outfile.close();
}

double func(double x){
	return x-1;
}