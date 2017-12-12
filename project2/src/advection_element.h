#ifndef ADVECTION_ELEMENT
#define ADVECTION_ELEMENT

#include <vector>

using namespace std;

/*
AdvectionElement handles data storage for the elements
of a finite element mesh. Elements are connected by 
neighbours through their pointers, and contain data on
the position on the domain as well as the value of the
solution on the boundary elements.
*/
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

	// Calculate the flux for advection equation
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

	// Numerical flux function for advection equation.
	// Uses LLF flux to approximate the flux numerically.
	virtual double h(double a, double b) {
		return flux(a);
	}
};

/*
BurgerElement inherits from AdvectionElement.
Both classes have identical functionality, but
BurgerElement overrides the flux and numerical
flux functions for its own calculations.
*/
class BurgerElement : public AdvectionElement {
public:
	// Flux function for Burgers' equation.
	double flux(double u){
		return 0.5*u*u;
	}

	// Numerical flux function for Burgers' equation.
	// Uses LLF flux to approximate the flux numerically.
	double h(double a, double b){
		return 0.5*(flux(a) + flux(b)) - 0.5*abs(std::max(a, b))*(b - a);
	}
};

#endif