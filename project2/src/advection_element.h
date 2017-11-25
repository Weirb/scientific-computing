#ifndef ADVECTION_ELEMENT
#define ADVECTION_ELEMENT

#include <vector>

using namespace std;

class AdvectionElement {
public:
	// Pointer to the left neighbour
	AdvectionElement *Left_neighbour_pt;

	// Pointer to the right neighbour
	AdvectionElement *Right_neighbour_pt;

	// Storage for the coordinates
	vector<double> X;

	// Storage for the unknowns
	vector<double> U;

	// Constructor: initialise the vectors to hold two entries.
	AdvectionElement() {
		X = vector<double>(2);
		U = vector<double>(2);
	}

	// Return the value of the coordinate at local coordinate s using
	// equation (1.2)
	double interpolated_x(double s) {
		return X[0] * 0.5*(1 - s) + X[1] * 0.5*(1 + s);
	}

	// Return the value of the unknown at local coordinate s using
	// equation (1.4)
	double interpolated_u(double s) {
		return U[0] * 0.5*(1 - s) + U[1] * 0.5*(1 + s);
	}
};

#endif