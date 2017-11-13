#include <iostream>
#include "mvector.h"
#include "mmatrix.h"
#include "mbandedmatrix.h"
#include "functions.h"

using namespace std;

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
