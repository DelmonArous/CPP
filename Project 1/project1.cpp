#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "time.h"
#include "armadillo"

using namespace std;
using namespace arma;

void initialize(int &n);
void initialize_vectors(int n, double *x, double *h_step, double *a, double *b, double *c, double *b_tilde, vec &b_tilde_lu);
void fw_bw_subst(int n, double *a, double *b, double *c, double *v, double *b_tilde);
void output(int n, double *x, double *h_step, double *v);

ofstream outfile1, outfile2;

int main() {

    // Declare variables
    int n;
    double *x, *h_step, *a, *b, *c, *v, *b_tilde;
	clock_t start1, start2, finish1, finish2; // declare start and final time

    // Read in output file, abort if there are too few command-line arguments
    outfile1.open("data_plot.txt", ios::out);
    outfile2.open("data.txt", ios::out);

    // Read in input data from screen
    initialize(n);
	
    // Allocate space in memory for the one-dimensional arrays h_step and computed_derivative
    x = new double[n];
    h_step = new double[n];
    a = new double[n];
    b = new double[n];
    c = new double[n];
    v = new double[n];
    b_tilde = new double[n];
    vec b_tilde_lu(n);

    // Initialize vectors needed
    initialize_vectors(n, x, h_step, a, b, c, b_tilde, b_tilde_lu);

    // Commence forward substitution, and then backward substitution
    start1 = clock();
	fw_bw_subst(n, a, b, c, v, b_tilde);
	finish1 = clock();

	cout << "Number of FLOPS for selfmade tridiagonal solver: " << 9*n << endl;
	cout << "Execution time for tridiagonal solver: " << ( (finish1 - start1)/float(CLOCKS_PER_SEC) ) << " s" <<endl;

    // LU decomposition initialization
	mat A(n,n), L(n,n), U(n,n);
	vec v_lu(n), y(n), off_diag(n-1);

	// Defining the tridiagonal matrix
	off_diag.fill(-1); // off-diagonal elements are defined in the vector off_diag of length (n-1)
	A = 2.0*eye(n,n); // diagonal elements of the tridiagonal matrix A are defines
	A.diag(1)= off_diag; A.diag(-1)= off_diag; // elements of off_diag are placed in the sub- and upper diagonal in A
	A(0,0) = 1.0; A(n-1,n-1) = 1.0; // taking boundary conditions into account by "expanding" A
	A(0,1) = 0.0; A(1,0)= 0.0; A(n-1,n-2)= 0.0; A(n-2,n-1) = 0.0;

	// LU decomposition
    start2 = clock();
	lu(L, U, A);
	y = inv(L)*b_tilde_lu;
	v_lu = inv(U)*y;
	finish2 = clock();

	cout << "Number of FLOPS for selfmade tridiagonal solver: " << floor(2.0/3.0*pow(n,3)) << endl;
	cout << "Execution time for LU decomposition: " << ( (finish2 - start2)/float(CLOCKS_PER_SEC)) << " s" << endl;

    // Print out results to file
    output(n, x, h_step, v);

    // Clear memory
    delete [] x;
    delete [] h_step;
    delete [] a;
    delete [] b;
    delete [] c;
    delete [] v;
    delete [] b_tilde;

    outfile1.close();
    outfile2.close();
    return 0;
}

void initialize(int &n) {
    cout << "Read in from screen:" << endl;
    cout << "\n Number of grid points: n = "; cin >> n;
    return;
}

void initialize_vectors(int n, double *x, double *h_step, double *a, double *b, double *c, double *b_tilde, vec &b_tilde_lu) {
    double h = 1.0 / (n-1.0);
    for(int i = 0; i < n; i++) {
        h_step[i] = h;
        x[i] = i*h;
        a[i] = -1;
        b[i] = 2;
        c[i] = -1;
        b_tilde[i] = h*h*100*exp(-10*x[i]);
		b_tilde_lu(i) = b_tilde[i];
    }
	b_tilde[0]=0.0; b_tilde[n-1]=0.0;
	b_tilde_lu(0)=0.0; b_tilde_lu(n-1)=0.0;
    a[0] = 0; c[n-1] = 0;
    return;
}

void fw_bw_subst(int n, double *a, double *b, double *c, double *v, double *b_tilde) {
    double m;
	// Forward substitution
    for(int i = 1; i < n-1; i++) { 
            m = a[i]/b[i-1]; // line a
            b[i] -= m*c[i-1]; // line b
            b_tilde[i] -= m*b_tilde[i-1]; // line c
        }
    // Backward substitution
    v[n-1] = b_tilde[n-1]/b[n-1]; // line d
    for(int i = n-2; i > 0; i--) {
        v[i] = (b_tilde[i] - c[i]*v[i+1])/b[i]; // line e
    }
    v[0] = 0; v[n-1] = 0;
    return;
}

void output(int n, double *x, double *h_step, double *v) {
    double u, x_temp, temp_error, rel_error = 0;
    outfile1 << setiosflags(ios::showpoint | ios::uppercase);
    outfile2 << setiosflags(ios::showpoint | ios::uppercase);

    for(int i = 1; i < n; i++) {
        u = 1.0-(1.0-exp(-10.0))*x[i]-exp(-10.0*x[i]);
        temp_error = log10(fabs( (v[i] - u)/u ));
        outfile1 << setw(15) << setprecision(8) << x[i];
        outfile1 << setw(15) << setprecision(8) << v[i];
        outfile1 << setw(15) << setprecision(8) << u << endl;
        outfile2 << setw(15) << setprecision(8) << log10(h_step[i]);
        outfile2 << setw(15) << setprecision(8) << temp_error << endl;
        if (fabs(temp_error) > fabs(rel_error)) {
            rel_error = temp_error; // update rel_error when temp_error is larger than rel_error
            x_temp = x[i]; // saves the position x_i for the update
        }
    }
    cout << "Maximum relative error: " << rel_error << " for x = " << x_temp << endl;
    return;
}