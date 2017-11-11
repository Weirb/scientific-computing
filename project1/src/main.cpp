#include <iostream>
#include "mvector.h"
#include "mmatrix.h"
#include "mbandedmatrix.h"
#include "functions.h"
#include <cmath>

using namespace std;

struct VectorSolution {
/*
VectorSolution struct explicitly defines the solution parameters.
These are the number of iterations taken to converge and the solution vector.
*/
	MVector solution;
	int iteration_count;
};

template<class Matrix>
VectorSolution cg(Matrix A, MVector b, MVector x0, int itercount, double tol){
/*
Conjugate Gradient
Compute the solution to A*x=b with initial guess x0.
Runs for no more than itercount iterations and finds solution to within tol.
*/
	int n = A.Rows();

	// Initial guess of the zero vector.
	MVector x = x0;

	MVector r = b - A*x;
	MVector p = r;
	MVector w, rk;
	double alpha = 0, beta = 0, rdotr = 0;

	int k = 0;
	for (; k < itercount; ++k){

		w = A * p;
		rdotr = dot(r, r);
		alpha =  rdotr / dot(p, w);

		// Update the current solution
		x = x + alpha*p;

		// Update the residual and check if norm is less than tolerance.
		// Exit loop if we have convergence.
		rk = r - alpha*w;
		if (rk.L2Norm() < tol){
			break;
		}

		// Compute the variable updates for p and r
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

int main(int argc, char** argv){

	int n = 5;
	MBandedMatrix A = laplacian_2d_banded(n);
	MVector b = create_vector1(n*n);
	VectorSolution sol = cg(A, b, MVector(n*n, 0), 1000, 1e-6);
	cout << sol.iteration_count << endl;
	cout << sol.solution << endl;

	MMatrix B = laplacian_2d(n);
	sol = cg(B, b, MVector(n*n, 0), 1000, 1e-6);
	cout << sol.iteration_count << endl;
	cout << sol.solution;

	return 0;
}
