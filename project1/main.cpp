#include <iostream>
#include "mvector.h"
#include "mmatrix.h"
#include "mbandedmatrix.h"
#include "functions.h"
#include <cmath>

using namespace std;

struct VectorSolution {
	MVector solution;
	int iteration_count;
};

template<class Matrix>
VectorSolution cg(Matrix A, MVector b, MVector x0, int itercount, double tol){

	int n = A.Rows();

	// Initial guess of the zero vector
	MVector x = x0;

	MVector s = b - A*x;
	MVector r = s;

	int k;
	for (k = 0; k < itercount; ++k){

		MVector w = A*s;

		double a = dot(r, s) / dot(w, s);

		x = x + a*s;

		MVector rk = r - a*w;

		if (r.L2Norm() < tol){
			break;
		}

		double p = pow(rk.L2Norm() / r.L2Norm(), 2);

		s = rk + p*s;

		r = rk;
	}

	return VectorSolution{x, k};
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
