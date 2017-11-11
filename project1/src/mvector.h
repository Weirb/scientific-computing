/*
File: mvector.h
MVector class represents a mathematical vector.
*/

#ifndef MVECTOR_H
#define MVECTOR_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using std::max;
using std::abs;
using std::vector;

class MVector {
public:
	// constructors
	MVector() {}
	MVector(int n) : v(n) {}
	MVector(int n, double x) : v(n, x) {}

	// access element (lvalue)
	double &operator[](int index) { return v[index]; }
	double &operator()(int index) { return v[index]; }

	// access element (rvalue)
	double operator[](int index) const { return v[index]; }
	double operator()(int index) const { return v[index]; }

	// number of elements
	int size() const { return v.size(); }

	// L2 norm, square root of the sum of squares of elements.
	double L2Norm() const {
		double norm = 0;
		for (int i = 0; i < size(); ++i)
			norm += pow(v[i], 2);

		return sqrt(norm);
	}

	// Linfinity norm, maximum absolute value of elements.
	double LInfNorm() const {
		double maxAbs = 0;
		for (int i = 0; i < size(); i++)
			maxAbs = max(abs(v[i]), maxAbs);

		return maxAbs;
	}

private:
	vector<double> v;
};

// Compute the dot product of two MVectors
inline double dot(const MVector& lhs, const MVector& rhs) {
	// Ensure operation is valid 
	if (lhs.size() != rhs.size()) {
		cerr << "Error in function: double dot(const MVector&, const MVector&);\n"
			 << "MVectors incompatible length.\n";
		exit(-1);
	}

	double t = 0;
	for (int i = 0; i < lhs.size(); ++i)
		t += lhs[i] * rhs[i];

	return t;
}

// Operator overload for "scalar * MVector"
inline MVector operator*(const double& lhs, const MVector& rhs) {
	MVector temp(rhs);
	for (int i = 0; i < temp.size(); i++)
		temp[i] *= lhs;

	return temp;
}

// Operator overload for "MVector * scalar"
inline MVector operator*(const MVector& lhs, const double& rhs) {
	// Switch order of operands to invoke operator*().
	return operator*(rhs, lhs);
}

// Operator overload for "MVector / scalar"
inline MVector operator/(const MVector& lhs, const double& rhs) {
	return operator*(lhs, 1 / rhs);
}

// Operator overload for "MVector + MVector"
inline MVector operator+(const MVector& lhs, const MVector& rhs) {
	// Ensure operation is valid 
	if (lhs.size() != rhs.size()) {
		cerr << "Error in function: MVector operator+(const MVector&, const MVector&);\n"
			 << "MVectors incompatible length.\n";
		exit(-1);
	}

	MVector temp(rhs);
	for (int i = 0; i < temp.size(); i++)
		temp[i] += lhs[i];

	return temp;
}

// Operator overload for "MVector - MVector"
inline MVector operator-(const MVector& lhs, const MVector& rhs) {
	// Ensure operation is valid 
	if (lhs.size() != rhs.size()) {
		cerr << "Error in function: MVector operator-(const MVector&, const MVector&);\n"
			 << "MVectors incompatible length.\n";
		exit(-1);
	}

	MVector temp(lhs);
	for (int i = 0; i < temp.size(); i++)
		temp[i] -= rhs[i];

	return temp;
}

// Overload operator for MVector to an output stream.
inline ostream& operator<<(ostream& stream, const MVector& v) {
	// Output each element on a new line.
	// No new line after the final element.
	int n = v.size();
	for (int i = 0; i < n - 1; ++i) {
		stream << v[i] << endl;
	} stream << v[n - 1];
	return stream;
}

#endif