#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <sstream>
#include "mvector.h"
#include "mmatrix.h"
#include "mbandedmatrix.h"
#include "advection_element.h"

using namespace std;

int get_element_count(const AdvectionElement* e) {
	AdvectionElement *current = e->Right_neighbour_pt;
	int c = 1;
	while (current != e) {
		c++;
		current = current->Right_neighbour_pt;
	}
	return c;
}

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

AdvectionElement* create_burger_elements(int N, double x0, double x1, double(*u)(double x)) {

	double h = (x1 - x0) / N;

	// Advection element is basically a circular doubly linked list
	AdvectionElement* start = new BurgerElement;
	start->X[0] = x0;
	start->X[1] = x0 + h;
	start->U[0] = u(start->X[0]);
	start->U[1] = u(start->X[1]);
	start->Left_neighbour_pt = start->Right_neighbour_pt = start;

	for (int i = 1; i < N; ++i) {
		// Keep track of the last element in the list
		AdvectionElement* last = start->Left_neighbour_pt;

		// Create the new advection element
		AdvectionElement* e = new BurgerElement;
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

void create_data_file() {
	double x0 = 0, x1 = 2 * M_PI;
	auto u = [](double x) { return 1.5 + sin(x); };

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

			f << x->interpolated_x(0) << "," << x->interpolated_u(0) << endl;
			x = x->Right_neighbour_pt;
		}

		f.close();
		f.clear();
	}
}

void create_animation_data(AdvectionElement* e) {

	ofstream f;
	f.open("animation_data.dat");
	if (!f.is_open()) {
		cout << "Could not open file." << endl;
		exit(-1);
	}

	int n = get_element_count(e);

	for (int k = 0; k < 48; ++k){

		e->timestep(0.05);

		for (int i = 0; i < n; ++i) {

			f << e->interpolated_u(0);
			if (i < n - 1)
				f << ",";
			e = e->Right_neighbour_pt;
		} f << endl;

	}

	cout << "Finished." << endl;

	f.close();
}

void create_data_file2(AdvectionElement* e) {

	ofstream f;
	f.open("plot_data.dat");
	if (!f.is_open()) {
		cout << "Could not open file." << endl;
		exit(-1);
	}

	int n = get_element_count(e);

	e->timestep(0.1);
	e->timestep(0.1);

	for (int i = 0; i < n; ++i) {

		f << e->interpolated_x(0) << "," << e->interpolated_u(0) << endl;
		e = e->Right_neighbour_pt;
	}

	f.close();
}

int main(int argc, char** argv){
	auto square_ic = [](double x) { return (0 <= x) && (x <= 1) ? 1. : 0.; };
	auto sine_ic = [](double x) { return 1.5 + sin(x); };

	AdvectionElement *b = create_advection_elements(50, 0, 2 * M_PI, square_ic);

	create_animation_data(b);

	// b->timestep(0.05);
	// b->timestep(0.05);
	// b->timestep(0.05);

	// AdvectionElement* c = b;
	// do {
	// 	cout << c->interpolated_x(0) << ", " << c->interpolated_u(0) << endl;
	// 	c = c->Right_neighbour_pt;
	// } while (b != c);

	// cin.get();
	return 0;
}