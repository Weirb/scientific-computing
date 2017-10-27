/*
Class that represents a mathematical matrix
*/

#ifndef MMATRIX_H
#define MMATRIX_H

#include <vector>
#include <iostream>

class MMatrix {
public:
	// constructors
	MMatrix() : nRows(0), nCols(0) { }
	MMatrix(int n, int m, double x = 0) : nRows(n), nCols(m), A(n * m, x) { }

	// set all matrix entries equal to a double
	MMatrix &operator=(double x) {
		for (unsigned int i = 0; i < nRows * nCols; i++)
			A[i] = x;

		return *this;
	}

	// access element, indexed by (row, column) [rvalue]
	double operator()(int i, int j) const {
		return A[j + i * nCols];
	}

	// access element, indexed by (row, column) [lvalue]
	double &operator()(int i, int j) {
		return A[j + i * nCols];
	}

	// size of matrix
	int Rows() const { return nRows; }
	int Cols() const { return nCols; }

private:
	unsigned int nRows, nCols;
	std::vector<double> A;
};



inline MVector operator*(const MMatrix& A, const MVector& x){
	if (A.Cols() != x.size())
		exit(-1);

	int n = A.Rows();
	MVector v(n, 0.);

	for (int i = 0; i < n; ++i){
		for (int j = 0; j < A.Cols(); ++j){
			v[i] += A(i, j)*x(j);
		}
	}
	return v;
}

inline MMatrix operator+(const MMatrix& A, const MMatrix& B) {
	if (!(A.Cols() == B.Cols() && A.Rows() == B.Rows()))
		exit(-1);

	MMatrix C = A;

	for (int i = 0; i < C.Rows(); ++i) {
		for (int j = 0; j < C.Cols(); ++j) {
			C(i, j) += B(i, j);
		}
	}
	return C;
}

inline ostream& operator<<(ostream& stream, const MMatrix& A) {
	for (int i = 0; i < A.Rows(); ++i) {
		for (int j = 0; j < A.Cols(); j++) {
			stream.width(8);
			stream << A(i, j);
		} stream << endl;
	}
	return stream;
}

#endif