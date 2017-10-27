/*
Class that represents a mathematical vector
*/


#ifndef MVECTOR_H
#define MVECTOR_H

#include <iostream>
#include <vector>
#include <algorithm>

using std::cerr;

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

	double L2Norm() const {
		double norm = 0;
		for (int i = 0; i < size(); ++i)
			norm += pow(v[i], 2);

		return sqrt(norm);
	}

	double LInfNorm() const {
		double maxAbs = 0;
		for (int i = 0; i < size(); i++)
			maxAbs = std::max(std::abs(v[i]), maxAbs);

		return maxAbs;
	}

private:
	std::vector<double> v;

};

double dot(const MVector& lhs, const MVector& rhs) {
	// Ensure operation is valid 
	if (lhs.size() != rhs.size()) {
		cerr << "Error in function: dot();\nVectors incompatible length.\n";
		exit(-1);
	}

	double t = 0;
	for (int i = 0; i < lhs.size(); ++i)
		t += lhs[i] * rhs[i];

	return t;
}

// Operator overload for "scalar * vector"
inline MVector operator*(const double& lhs, const MVector& rhs) {
	MVector temp(rhs);
	for (int i = 0; i < temp.size(); i++)
		temp[i] *= lhs;

	return temp;
}

// Operator overload for "vector * scalar"
inline MVector operator*(const MVector& lhs, const double& rhs) {
	MVector temp(lhs);
	for (int i = 0; i < temp.size(); i++)
		temp[i] *= rhs;

	return temp;
}

// Operator overload for "vector / scalar"
inline MVector operator/(const MVector& lhs, const double& rhs) {
	return operator*(lhs, 1 / rhs);
}

// Operator overload for "vector + vector"
inline MVector operator+(const MVector& lhs, const MVector& rhs) {
	// Ensure operation is valid 
	if (lhs.size() != rhs.size()) {
		cerr << "Error in function: MVector operator+();\nVectors incompatible length.\n";
		exit(-1);
	}

	MVector temp(rhs);
	for (int i = 0; i < temp.size(); i++)
		temp[i] += lhs[i];

	return temp;
}

// Operator overload for "vector - vector"
inline MVector operator-(const MVector& lhs, const MVector& rhs) {
	// Ensure operation is valid 
	if (lhs.size() != rhs.size()) {
		cerr << "Error in function: MVector operator-();\nVectors incompatible length.\n";
		exit(-1);
	}

	MVector temp(lhs);
	for (int i = 0; i < temp.size(); i++)
		temp[i] -= rhs[i];

	return temp;
}

inline ostream& operator<<(ostream& stream, const MVector& v) {
	for (int i = 0; i < v.size() - 1; ++i) {
		stream << v[i] << endl;
	} stream << v[v.size() - 1];
	return stream;
}

#endif