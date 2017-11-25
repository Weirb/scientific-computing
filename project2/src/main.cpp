#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <sstream>
#include "mvector.h"
#include "mmatrix.h"
#include "mbandedmatrix.h"
#include "advection_element.h"

using namespace std;

AdvectionElement* create_advection_elements(int N, double x0, double x1, double (*u)(double x)) {

	double h = (x1 - x0) / N;

	// Advection element is basically a circular doubly linked list
	AdvectionElement* start = new AdvectionElement;
	start->X[0] = x0;
	start->X[1] = x0 + h;
	start->U[0] = u(start->X[0]);
	start->U[1] = u(start->X[1]);
	start->Left_neighbour_pt = start->Right_neighbour_pt = start;

	for (int i = 1; i < N; ++i) {
		// Keep track of the last element in the list
		AdvectionElement* last = start->Left_neighbour_pt;

		// Create the new advection element
		AdvectionElement* e = new AdvectionElement;
		e->X[0] = x0 + i*h;
		e->X[1] = x0 + (i + 1)*h;
		e->U[0] = u(e->X[0]);
		e->U[1] = u(e->X[1]);

		// Update the pointers to maintain the correct sequence of elements
		// - e is the new last element, so left neighbour is last and right neighbour is start
		// - start's new left neighbour is e
		// - last's new right neighbour is e
		e->Right_neighbour_pt = start;
		e->Left_neighbour_pt = last;
		start->Left_neighbour_pt = e;
		last->Right_neighbour_pt = e;
	}

	return start;
}

int main(int argc, char** argv){

	double x0 = 0, x1 = 2 * M_PI;
	auto u = [](double x){ return 1.5 + sin(x); };
	
	vector<int> ns = { 10, 100, 200 };
	for (auto it = ns.begin(); it != ns.end(); ++it) {

		ofstream f;
		stringstream ss;
		ss << "plot_data_" << *it << ".dat";
		f.open(ss.str());
		if (!f.is_open()) {
			cout << "Could not open file." << endl;
			break;
		}

		AdvectionElement* x = create_advection_elements(*it, x0, x1, u);
		for (int i = 0; i < *it; ++i) {

			f << x->interpolated_x(0.5) << "," << x->interpolated_u(0.5) << endl;
			x = x->Right_neighbour_pt;
		}

		f.close();
		f.clear();
	}

	cin.get();
	return 0;
}