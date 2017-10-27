#include <iostream>
#include "mvector.h"
#include "mmatrix.h"
#include "mbandedmatrix.h"
#include "functions.h"
#include <cmath>

using namespace std;

typedef pair < MVector, int > vectorpair;



MMatrix kronecker(MMatrix A, MMatrix B) {

	int m = A.Rows(); int n = A.Cols();
	int p = B.Rows(); int q = B.Cols();

	MMatrix C(m*p, n*q);

	for (int i = 0; i < m*p; ++i) {
		for (int j = 0; j < n*q; ++j) {
			C(i, j) = A(floor(i / p), floor(j / q)) * B(i - floor(i / p)*p, j - floor(j / q)*q);
		}
	}

	return C;
}




vectorpair cg(MMatrix A, MVector b, MVector x0, int itercount, double tol){

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

	return vectorpair(x, k);
}

vectorpair cg(MBandedMatrix A, MVector b, MVector x0, int itercount, double tol) {

	int n = A.Rows();

	// Initial guess of the zero vector
	MVector x = x0;

	MVector s = b - A*x;
	MVector r = s;

	int k;
	for (k = 0; k < itercount; ++k) {

		MVector w = A*s;

		double a = dot(r, s) / dot(w, s);

		x = x + a*s;

		MVector rk = r - a*w;

		if (r.L2Norm() < tol) {
			break;
		}

		double p = pow(rk.L2Norm() / r.L2Norm(), 2);

		s = rk + p*s;

		r = rk;
	}

	return vectorpair(x, k);
}

//void cg_test() {
//
//	int n = 25;
//	double tol = 1e-6;
//
//	MMatrix A = create_matrix1(n);
//	MVector b = create_vector2(n);
//	MVector x0(n, 0);
//	vectorpair v = cg(A, b, x0, 1000, tol);
//	MVector x = v.first; int count = v.second;
//
//	cout << A << endl;
//
//	MVector r = A*x - b;
//	cout << r.L2Norm() << endl;
//	cout << count;
//
//}

int main(int argc, char** argv){

	// Kronecker test
	//MMatrix A(2, 2); 
	//MMatrix B(2, 2);
	//A(0, 0) = 1; A(0, 1) = 2;
	//A(1, 0) = 3; A(1, 1) = 4;
	//B(0, 0) = 0; B(0, 1) = 5;
	//B(1, 0) = 6; B(1, 1) = 7;
	//cout << kronecker(A, B);

	int n = 100;
	MBandedMatrix A = laplacian_2d_banded(n);
	MVector b = create_vector1(n*n);
	vectorpair x = cg(A, b, MVector(n*n, 0), 1000, 1e-6);
	cout << x.first << endl;
	cout << x.second;

	cin.get();

	return 0;




	//int n = 4;
	//MBandedMatrix A(n,n, 1, 1);
	//MMatrix M = create_matrix1(n);
	//for (int i = 0; i < n; ++i) {
	//	int jmin = max(min(i - A.LBands(), A.Cols()), 0);
	//	int jmax = min(i + A.RBands() + 1, A.Cols());

	//	for (int j = jmin; j < jmax; j++) {
	//		A(i, j) = M(i, j);
	//	}
	//}
	//cout << A << endl;
	

	//int n = pow(2,12);
	//double tol = 1e-6;
	//MBandedMatrix A = create_matrix1(n);
	//MVector b = create_vector1(n);

	//MVector x0(n, 0);
	//vectorpair v = cg(A, b, x0, 1000, tol);
	//MVector x = v.first; int count = v.second;

	//MVector r = A*x - b;
	//cout << r.L2Norm() << endl;
	//cout << count;



	//cout << A << endl;
	//cout << x << endl << endl;
	//cout << A*x;

	cin.get();
	return 0;
}



//void test_operators() {
//	// TEST OPERATORS
//	double v0[] = { 0.1, 4.8, 3.7 };
//	double w0[] = { 3.1, 8.5, 3.6 };
//	double x0[] = { 5.8, 7.4, 12.4 };
//
//	MVector v(3, v0);
//	MVector w(3, w0);
//	MVector x(3, x0);
//
//	MVector u = 4.7*v + 1.3*w - 6.7*x;
//
//	cout << u << endl;
//}
//
//void test_norms(){
//	// TEST NORMS, DOT PRODUCT
//	double u0[] = { 1.5, 1.3, 2.8 };
//	double v0[] = { 6.5, 2.7, 2.9 };
//	double w0[] = { 0.1, -7.2, 3.4 };
//
//	MVector u(3, u0);
//	MVector v(3, v0);
//	MVector w(3, w0);
//
//
//	if (u.L2Norm() - 3.4322 < 1e-6)
//		cout << "u correct" << endl;
//	if (v.L2Norm() - 7.61249 < 1e-6)
//		cout << "v correct" << endl;
//	if (w.L2Norm() - 7.96304 < 1e-6)
//		cout << "w correct" << endl;
//
//	if (dot(u, u) / dot(v, w) - (-1.31915) < 1e-4)
//		cout << "dot correct" << endl;
//}

void test_matrix1() {

	MMatrix m(4, 3);

	m(0, 0) = 1;
	m(0, 1) = 2;
	m(0, 2) = 3;

	m(1, 0) = 2;
	m(1, 1) = 3;
	m(1, 2) = 4;

	m(2, 0) = 3;
	m(2, 1) = 4;
	m(2, 2) = 5;

	m(3, 0) = 4;
	m(3, 1) = 5;
	m(3, 2) = 6;

}

void test_matrix2() {

	/* MATLAB TEST
	m = 4; n = 3;
	A = zeros(m, n);
	x = [0.5, 1.6, 3.2]';

	for i = 1:m
	for j = 1:n
	A(i,j) = 3*(i-1) + (j-1);
	end
	end

	A*x
	*/

	int m = 4, n = 3;
	MMatrix A(m, n);

	for (int i = 0; i < m; ++i){
		for (int j = 0; j < n; ++j){
			A(i, j) = 3 * i + j;
		}
	}

	MVector x(3);
	x(0) = 0.5;
	x(1) = 1.6;
	x(2) = 3.2;

	MVector b = A*x;

	cout << b << endl;
}
