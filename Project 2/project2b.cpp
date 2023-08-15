#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "armadillo"
#include "time.h"

using namespace std;
using namespace arma;

// Function declaration
void initialize(int &n, double &rho_max, double &omega_r);
double potential1(double x);
double potential2(double x, double omega_r);
void initialize_arrays(int n, double rho_max, double omega_r, vec &d, vec &e, vec &rho, mat &A);
void jacobi_method(int n, mat &A, mat &R);
double max_offdiag(int n, mat &A, int &k, int &l);
void rotate(int n, mat &A, mat &R, int k, int l);
void output(int n, double rho_max, mat R, vec rho, vec lambda, vec eigval, clock_t start1, clock_t start2, clock_t finish1, clock_t finish2, mat prob_vectors);

ofstream outfile; // output file as global variable 
ofstream outfile_plot;

int main() {
	int n;
	double rho_max, omega_r;
	clock_t start1, start2, finish1, finish2; // declare start and finish time
	
	outfile.open("data_norep_w050.txt", ios::out);
	outfile_plot.open("data_plot_norep_w050.txt", ios::out);

	// Read in data
	initialize(n, rho_max, omega_r);
	
	// Initialize constants and arrays
	vec d(n), e(n-1), eigval(n), lambda(n), rho(n);
	mat A(n,n), R(n,n), eigvec(n,n), norm_R(n,n), prob_vectors(n,3);
	initialize_arrays(n, rho_max, omega_r, d, e, rho, A);
	
	// Solve the eigensystem using Armadillo function and store eigenvalues 
	start1 = clock();
	eig_sym(eigval, eigvec, A);
	finish1 = clock();
	
	eigval = sort(eigval); // sort out the elements in the vector eigval

	// Jacobi method
	start2 = clock();
	jacobi_method(n, A, R);
	finish2 = clock();
	
	// Sort out and store the eigenvalues in the vector lambda
	lambda = sort(A.diag(0));

	//finding the lowest eigenvalues
	uvec col_indices;
	col_indices= sort_index(A.diag(0));

	// Putting the three lowest eigenvector in the matrix prob_vectors
	// Have used fixed end points so two of the eigenvalues are 1 so, we dont use them
	for(int i=0; i<n; i++){
		for(int j=0; j<=2; j++){
			prob_vectors(i,j) = R(i,col_indices(j+2))*R(i,col_indices(j+2));
		}
	}

	// We normalize the lowest eigenvectors before plotting
	vec norm_factor(3);
	norm_factor.fill(0.0);
	for(int j=0; j<=2; j++){
		for(int i=0; i<n; i++){
			norm_factor(j) += prob_vectors(i,j)*prob_vectors(i,j);
		}
		norm_factor(j) = sqrt(norm_factor(j));
	}

	// Here we divide by the norm to actually normalize the vectors
	for(int i=0; i<n; i++){
		for(int j=0; j<=2; j++){
			prob_vectors(i,j) = prob_vectors(i,j)/norm_factor(j);
		}
	}

	// Write to output file
	output(n, rho_max, R, rho, lambda, eigval, start1, start2, finish1, finish2, prob_vectors);
	outfile_plot.close();
	outfile.close();
	return 0;
}

void initialize(int &n, double &rho_max, double &omega_r) {
	// Printing out to screen what values need to be set for the program to 
	// run and reading these values
	cout << "Enter maximum value of rho: ";
	cin >> rho_max;
	cout << "\nEnter number of grid points: n = ";
	cin >> n;
	cout << "\nEnter frequency of potential: ";
	cin >> omega_r;
	return;
}

double potential1(double x) {
	// Harmonic oscillator potential
	return x*x;
}
double potential2(double x, double omega_r) {
	// Harmonic oscillator potential with additional term for electron repulsion
	return omega_r*omega_r*x*x + (1.0/x);
}
double potential3(double x, double omega_r) {
	// Harmonic oscillator potential without additional term for electron repulsion
	return omega_r*omega_r*x*x;
}

void initialize_arrays(int n, double rho_max, double omega_r, vec &d, vec &e, vec &rho, mat &A) {
	A.zeros();
	double h = rho_max/n;
	e.fill( (-1.0/(h*h)) );

	// Do not include the boundaries, these are set by default
	for(int i = 1; i < n-1; i++) { 
		rho(i) = i*h;
		d(i) = (2.0/(h*h)) + potential3(rho(i),omega_r); // or + potential2(rho(i), omega_r) depending on the eigenproblem
	}
	rho(0) = 0.0; rho(n-1) = rho_max;
	d(0)= 1.0; d(n-1)= 1.0; // boundary conditions for the radial wavefunction yields u(0)=u(rho_max)=0
	e(0)= 0.0; e(n-2)= 0.0; // therefore the vectors d,e inherit the same values at the boundaries
	A.diag(0) = d; A.diag(-1) = e; A.diag(1) = e;

	return;
}

void jacobi_method(int n, mat &A, mat &R) {
	R.eye(); // setting up the eigenvector matrix, this will actually give us a normalized matrix afer all the rotations
	int row_max, col_max, i = 0;

	// setting the precision criteria and max number of iterations we wanna make
	double epsilon = 1.0E-9;
	double A_max = max_offdiag(n, A, row_max, col_max);
	double max_number_iterations = (double) n * (double) n * (double) n;
	while( A_max*A_max > epsilon && (double) i < max_number_iterations ) {
		A_max = max_offdiag(n, A, row_max, col_max);
		rotate(n, A, R, row_max, col_max); // S^T * A * S
		i++;
	}
	outfile << "Number of iterations: " << i << "\n";
	
	return; 
}

double max_offdiag(int n, mat &A, int &k, int &l) {
	// setting max value to 0 at the begining so we dont get something random from the memory
	double max = 0.0;

	// Double for-loop to run over every matrix element (this can be optimized for tri-diagonal matrix)
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			// We only want the off-diagonal to be 0 therefore we do not include the diagonal
			if ( (i != j) && (fabs(A(i,j)) > max) ) {
				k = i; l = j;
				max = fabs(A(i,j));
			}
		}
	}
	return max;
}

void rotate(int n, mat &A, mat &R, int k, int l) {
	double s, c;

	// If A(k,l) = 0.0 we shouldn't be i this loop as this is lower than our precision
	if ( A(k,l) != 0.0 ) {
		double t, tau;
		tau = ( A(l,l) - A(k,k))/(2*A(k,l));

		// choosing the lowest of the roots
		if ( tau > 0 ) { 
			t = -1.0*tau + sqrt(1.0 + tau*tau);
		}
		else {
			t = -1.0*tau - sqrt(1.0 + tau*tau);
		}
		c = 1.0/sqrt(1.0+t*t);
		s = t*c;
	} 
	// In case the while loop in the Jacobis method does not work we hard-code this in
	else {
		c = 1.0;
		s = 0.0;
	}

	double a_kk, a_kl, a_ll, a_ik, a_il, r_ik, r_il;
	a_kk = A(k,k);
	a_kl = A(k,l); 
	a_ll = A(l,l);

	// Changing the matrix elements with indices k and l
	A(k,k) = a_kk*c*c - 2.0*a_kl*c*s + a_ll*s*s;
	A(l,l) = a_ll*c*c + 2.0*a_kl*c*s + a_kk*s*s;
	A(k,l) = 0.0; A(l,k) = 0.0; // hard-coding of the zeros to save flops

	// Changing the remaining elements
	for ( int i = 0; i < n; i++ ) {
		if ( i != k && i != l ) {
			a_ik = A(i,k); a_il = A(i,l);
			A(i,k) = a_ik*c - a_il*s;
			A(k,i) = A(i,k);
			A(i,l) = a_il*c + a_ik*s;
			A(l,i) = A(i,l);
		}
		// Compute the new eigenvectors
		r_ik = R(i,k); r_il = R(i,l);
		R(i,k) = r_ik*c - r_il*s;
		R(i,l) = r_il*c + r_ik*s;
	}
	return;
}

void output(int n, double rho_max, mat R, vec rho, vec lambda, vec eigval, clock_t start1, clock_t start2, clock_t finish1, clock_t finish2,mat prob_vectors) {
	vec exact_eigenval(3); 
	exact_eigenval(0) = 3.0; exact_eigenval(1) = 7.0; 	exact_eigenval(2) = 11.0; 
	outfile << setiosflags(ios::showpoint | ios::uppercase);
	outfile << "Rho_max = " << setw(15) << setprecision(8) << rho_max << endl;
	outfile << "Number of steps: n = " << setw(15) << n << endl;
	outfile << "Execution time by using Armadillo function: " << ( float(finish1 - start1) / CLOCKS_PER_SEC ) << " s" << endl; 
	outfile << "Execution time by using Jacobis method: " << ( float(finish2 - start2) / CLOCKS_PER_SEC ) << " s" << endl;
	outfile << "Three lowest eigenvalues: " << endl;
	for(int i = 2; i < 5; i++) {
		outfile << setw(15) << setprecision(8) << lambda(i) << "\t";
		outfile << setw(15) << setprecision(8) << eigval(i) << endl;
	}
	outfile << "The relative errors: " << endl;
	for(int i = 2; i < 5; i++) {
		outfile << setw(15) << setprecision(8) << log10( fabs((lambda(i)-exact_eigenval(i-2))/exact_eigenval(i-2)) ) << "\t";
		outfile << setw(15) << setprecision(8) << log10( fabs((eigval(i)-exact_eigenval(i-2))/exact_eigenval(i-2)) ) << endl;
	}
	for(int i = 0; i<n; i++){
		outfile_plot<<rho(i)<<"\t";
		for(int j = 0; j<3; j++){
			outfile_plot<< setw(15) <<prob_vectors(i,j)<<"\t";
		}
		outfile_plot<<endl;
	}
	return;
	
}