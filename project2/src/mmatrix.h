/*
File: mmatrix.h
MMatrix class represents a mathematical matrix.
*/

#ifndef MMATRIX_H
#define MMATRIX_H

#include <vector>
#include <iostream>

using namespace std;

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

	// size of matrix: rows and columns
	int Rows() const { return nRows; }
	int Cols() const { return nCols; }

private:
	unsigned int nRows, nCols;
	vector<double> A;
};

// Overload operator for MMatrix * MVector
inline MVector operator*(const MMatrix& A, const MVector& x){
	// Ensure operation is valid 
	if (A.Cols() != x.size()){
		cerr << "Error in function: MVector operator+(const MMatrix&, const MVector&);\n"
			 << "Operands incompatible size.\n";
		exit(-1);
	}

	int n = A.Rows();
	MVector v(n, 0.);

	for (int i = 0; i < n; ++i){
		for (int j = 0; j < A.Cols(); ++j){
			v[i] += A(i, j)*x(j);
		}
	}
	return v;
}

// Overload operator for MMatrix + MMatrix
inline MMatrix operator+(const MMatrix& A, const MMatrix& B) {
	// Ensure operation is valid 
	if (!(A.Cols() == B.Cols() && A.Rows() == B.Rows())){
		cerr << "Error in function: MVector operator+(const MMatrix&, const MMatrix&);\n"
			 << "Operands incompatible size.\n";
		exit(-1);
	}

	// Loop over all rows, columns of the matrix.
	// Element wise add the matrices together.
	MMatrix C = A;
	for (int i = 0; i < C.Rows(); ++i) {
		for (int j = 0; j < C.Cols(); ++j) {
			C(i, j) += B(i, j);
		}
	}
	return C;
}

// Overload operator for MMatrix to an output stream.
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