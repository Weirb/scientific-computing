#include "functions.h"

MMatrix laplacian_1d(int n) {

	MMatrix A(n, n, 0.);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (i == j)
				A(i, j) = 2;
			else if (abs(i - j) == 1)
				A(i, j) = -1;
		}
	}

	return A;
}

MBandedMatrix laplacian_1d_banded(int n) {

	MBandedMatrix A(n, n, 1, 1, 0.);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (i == j)
				A(i, j) = 2;
			else if (abs(i - j) == 1)
				A(i, j) = -1;
		}
	}

	return A;
}

MMatrix create_matrix2(int n, double m) {

	MMatrix A(n, n, 0.);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (i == j)
				A(i, j) = 2 * pow(i + 1, 2) + m;
			else if (abs(i - j) == 1)
				A(i, j) = -pow(i + 1, 2);
		}
	}

	return A;
}

MMatrix laplacian_2d(int n) {

	MMatrix A(n*n, n*n, 0.);

	//MMatrix T = create_matrix1(n);
	//MMatrix I(n, n);
	//for (int i = 0; i < n; ++i)
	//	I(i, i) = 1;

	//A = kronecker(I, T) + kronecker(T, I);
	for (int i = 0; i < n*n; ++i) {
		for (int j = 0; j < n*n; ++j) {
			if (i == j)
				A(i, j) = 4;
			else if (abs(i - j) == n)
				A(i, j) = -1;
			else if (abs(i - j) == 1 && (i + j) % (2 * n) != 2 * n - 1)
				A(i, j) = -1;
		}
	}
	return A;
}

MBandedMatrix laplacian_2d_banded(int n) {

	MBandedMatrix A(n*n, n*n, n, n);

	for (int i = 0; i < n*n; ++i) {
		for (int j = 0; j < n*n; ++j) {
			if (i == j)
				A(i, j) = 4;
			else if (abs(i - j) == n)
				A(i, j) = -1;
			else if (abs(i - j) == 1 && (i + j) % (2 * n) != 2 * n - 1)
				A(i, j) = -1;
		}
	}

	return A;
}

MVector create_vector1(int n) {
	return MVector(n, 1 / pow(n + 1, 2));
}

MVector create_vector2(int n) {
	return MVector(n, 2.5);
}