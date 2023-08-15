#include <iostream>
#include <cstdlib>
#include "time.h"
#include "armadillo"

using namespace std;
using namespace arma;

int main() {

 int N = 2000; // dimension of matrices
 clock_t start1, start2, finish1, finish2; // declare start and final time

 // Memory allocation of matrices
 double **A, **B, **C;
 
 A = new double*[N];
 B = new double*[N];
 C = new double*[N];
 
 for(int i = 0; i < N; i++) {
	 A[i] = new double[N];
	 B[i] = new double[N];
	 C[i] = new double[N];
 }
 
 // Initialization of matrices
 for(int i = 0; i < N; i++) {
	 for(int j = 0; j < N; j++) {
		 A[i][j] = 0.0; // initialize all elements to zero
		 B[i][j] = ((double) rand() / (RAND_MAX + 1)); // initialize all elements to be random
		 C[i][j] = ((double) rand() / (RAND_MAX + 1)); // initialize all elements to be random
	 }
 }

 // Matrix-matrix multiplication
 start1 = clock();
 for(int k = 0; k < N; k++) {
	 for(int i = 0; i < N; i++) {
		 for(int j = 0; j < N; j++) {
			 A[i][j] += B[i][k]*C[k][j];
		 }
	 }
 }
 finish1 = clock();
 cout<<"Execution time for tridiagonal solver: " << ( (finish1 - start1)/CLOCKS_PER_SEC )<< " s" << endl;

// Free space by deleting the matrices
 for(int i = 0; i < N; i++) {
	 delete [] A[i];
	 delete [] B[i];
	 delete [] C[i];
 }
 delete [] A;
 delete [] B;
 delete [] C;

 // Declarations of variables 
 mat A_;
 mat B_ = randu(N,N);
 mat C_ = randu(N,N);

 start2 = clock();
 A_ = B_*C_;
 finish2 = clock(); 
 
 cout<<"Execution time for Library Solver: "<< ( (finish2 - start2)/CLOCKS_PER_SEC ) << " s" << endl;

 return 0;
}