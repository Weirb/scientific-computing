/*

*/

#ifndef MBANDEDMATRIX_H
#define MBANDEDMATRIX_H

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

class MBandedMatrix {
public:
	// constructors
	MBandedMatrix() : nRows(0), nCols(0) {}
	MBandedMatrix(int n, int m, int lband, int rband, double x = 0.) :
		nRows(n), nCols(m), A(n, lband + rband + 1, x), l(lband), r(rband) {}

	// access element [rvalue]
	double operator()(int i, int j) const {
		//cout << "Accessing: " << j + l - i + i * nCols << endl;
		//cout << "Size: " << A.size() << endl << endl;

		return A(i, j + l - i);
	}

	// access element [lvalue]
	double &operator()(int i, int j) {
		return A(i, j + l - i);
	}

	// size of matrix
	int Rows() const { return nRows; }
	int Cols() const { return nCols; }

	int Bands() const { return r + l + 1; } // total number of bands
	int LBands() const { return l; } // number of left bands
	int RBands() const { return r; } // number of right bands

private:
	MMatrix A;
	unsigned int nRows, nCols;
	int l, r; // number of left/right diagonals
};

inline MVector operator*(const MBandedMatrix A, const MVector& x) {

	int n = A.Rows();
	MVector y(n, 0.);

	for (int i = 0; i < n; ++i) {
		int jmin = max(min(i - A.LBands(), A.Cols()), 0);
		int jmax = min(i + A.RBands() + 1, A.Cols());

		for (int j = jmin; j < jmax; j++) {
			y(i) += A(i, j) * x(j);
		}
	}
	return y;
}

inline ostream& operator<<(ostream& output, const MBandedMatrix& banded) {
	int r = banded.Rows(), c = banded.Cols();
	for (int i = 0; i < banded.Rows(); i++) {
		// calculate position of lower and upper band
		int jmin = max(min(i - banded.LBands(), banded.Cols()), 0);
		int jmax = min(i + banded.RBands() + 1, banded.Cols());

		for (int j = 0; j < jmin; j++) {
			output.width(8);
			output << 0;
		}
		for (int j = jmin; j < jmax; j++) {
			output.width(8);
			output << banded(i, j);
		}
		for (int j = jmax; j < c; j++) {
			output.width(8);
			output << 0;
		}
		std::cout << endl;
	}
	return output;
}

#endif
