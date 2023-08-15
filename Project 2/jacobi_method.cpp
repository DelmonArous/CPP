#include <iostream>
#include <cmath>
#include "armadillo"

using namespace std;
using namespace arma;

void jacobi_method(int n, mat &A, mat &R) {
	
	// Setting up the eigenvector matrix
	R.eye();

	int k, l, i=0;
	uword row_max, col_max; 
	double A_max = fabs(A.max(row_max, col_max));
	double epsilon = 1.0E-7;
	double max_number_iterations = (double) n * (double) n * (double) n;
	
	while( fabs(A_max) > epsilon && (double) i < max_number_iterations ) {
		A_max = fabs(A.max(row_max, col_max));
		rotate(A, R, row_max, col_max, n);
		i++;
	}
	cout << "Number of iterations: " << i << "\n";
	
	return; 
}



