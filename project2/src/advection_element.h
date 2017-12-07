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

	// Temporary storage for updating the unknown
	MVector U_update;

	// Constructor: initialise the vectors to hold two entries.
	AdvectionElement() {
		X = vector<double>(2);
		U = vector<double>(2);
		U_update = MVector(2, 0.);
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

	// Calculate the flux
	virtual double flux(double u) {
		return u;
	}

	// Calculate the integral of the flux function over the element
	// using the two-point Gauss rule
	double integrate_flux() {
		// Let's not use the sqrt function here
		// Approximate value of 1/sqrt(3) would be better to avoid
		// expensive function calls. Take 16 digits.
		double one_over_sqrt3 = 0.5773502691896257;
		return flux(interpolated_u(-one_over_sqrt3)) + flux(interpolated_u(one_over_sqrt3));
	}

	virtual double h(double a, double b) {
		// return 0.5 * (flux(a) + flux(b));
		return flux(a);
	}

	void timestep(double dt) {

		AdvectionElement* start = this;
		AdvectionElement* e = start;

		// Create the inverse of the mass matrix
		MMatrix Minv(2, 2, 0.);
		Minv(0, 0) = 2; Minv(0, 1) = -1;
		Minv(1, 0) = -1; Minv(1, 1) = 2;
		Minv = Minv * (2. / (e->X[1] - e->X[0]));

		// Compute the update values
		do {
			
			// Create temporary storage for the updated values
			e->U_update[0] = e->U[0];
			e->U_update[1] = e->U[1];

			// Flux vector
			MVector F(2);
			F(0) = -0.5; 
			F(1) = 0.5;
			F = F * e->integrate_flux();

			// Integral of the flux vector
			MVector H(2);
			H(0) = e->h(e->Left_neighbour_pt->U[1], e->U[0]);
			H(1) = -e->h(e->U[1], e->Right_neighbour_pt->U[0]);
			
			// Perform the update step
			e->U_update = e->U_update + dt*Minv*(F+H);

			e = e->Right_neighbour_pt;

		} while (e != start);

		// Update the values into U
		do {
			e->U[0] = e->U_update[0];
			e->U[1] = e->U_update[1];
			e = e->Right_neighbour_pt;
		} while (e != start);
	}
};

#endif