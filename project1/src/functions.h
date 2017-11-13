#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "mvector.h"
#include "mmatrix.h"
#include "mbandedmatrix.h"

struct VectorSolution {
/*
VectorSolution struct explicitly defines the solution parameters.
These are the number of iterations taken to converge and the solution vector.
*/
	MVector solution;
	int iteration_count;
};

// Templated function needs to be kept in header
template <class Matrix>
inline VectorSolution cg(const Matrix& A, 
				  const MVector& b, 
				  const MVector& x0, 
				  int max_iter_count, 
				  double tol) {
/*
Conjugate gradient algorithm
Compute the solution to A*x=b,
Inputs:
	Matrix A: 					Symmetric pos-def matrix, coefficients
	MVector b:					Right hand side of equation, data
	MVector x0: 				Initial guess
	int 	max_iter_count: 	Maximum number of iterations to compute
	double 	tol:				Error tolerance to find solution
Outputs:
	VectorSolution sol:			Contains solution and number of iterations 
								taken to converge to within tolerance
*/

	// Initialisation of variables
	int n = A.Rows();
	MVector x = x0; 		// Solution vector
	MVector r = b - A*x;	// Residual vector
	MVector p = r;			// Conjugate search direction
	MVector w;				// Store result of A*p for reduced computation
	MVector rk; 			// Updated residual vector
	double alpha = 0;		// Step size
	double beta = 0;		// Scalar to make updated p conjugate to the old p
	double rdotr = 0;		// Store result of dot(r,r) for reduced computation

	// Main loop, iteration counter k
	int k = 0;
	for (; k < max_iter_count; ++k){

		// Initialise reusable variables - avoid recomputations
		w = A * p;
		rdotr = dot(r, r);

		// Compute step size
		alpha =  rdotr / dot(p, w);

		// Update the current solution
		x = x + alpha*p;

		// Update the residual and check if norm is less than tolerance.
		// Exit loop if we have convergence.
		rk = r - alpha*w;
		if (rk.L2Norm() < tol){
			break;
		}

		// Compute the variable updates for beta, p, and r
		beta = dot(rk, rk) / rdotr;
		p = rk + beta * p;
		r = rk;
	}

	// Create the solution struct with variables to return.
	VectorSolution sol;
	sol.iteration_count = k;
	sol.solution = x;

	return sol;
}

MMatrix laplacian_1d(int n);
MBandedMatrix laplacian_1d_banded(int n);

MMatrix create_matrix2(int n, double m);

MMatrix laplacian_2d(int n);
MBandedMatrix laplacian_2d_banded(int n);

MVector create_vector1(int n);
MVector create_vector2(int n);

MMatrix kronecker(MMatrix A, MMatrix B);

#endif
