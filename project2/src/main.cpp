#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <sstream>
#include "mvector.h"
#include "mmatrix.h"
#include "mbandedmatrix.h"
#include "advection_element.h"
#include "mesh.h"

using namespace std;

template<class Element>
void create_data_file(Mesh<Element> &mesh, double dt, int N, string filename) {

	ofstream f(filename);
	if (!f.is_open()) {
		cout << "Could not open file." << endl;
		exit(-1);
	}

	for (int i = 0; i < N; ++i)
		mesh.timestep(dt);

	MeshData data = mesh.get_mesh_data(0);
	for (int i = 0; i < mesh.n_elements; ++i) {
		f << data.x[i] << "," << data.u[i]<< endl;
	}
		
	f.close();
}

template<class Element>
void create_animation_data(Mesh<Element> &mesh, double dt, int n_iterations, string filename) {

	ofstream f(filename);
	if (!f.is_open()) {
		cout << "Could not open file." << endl;
		exit(-1);
	}

	MeshData data;
	for (int k = 0; k < n_iterations; ++k){

		mesh.timestep(dt);
		data = mesh.get_mesh_data(0);
		for (int i = 0; i < mesh.n_elements; ++i) {
			f << data.u[i];
			if (i < mesh.n_elements - 1)
				f << ",";
		} f << endl;
	}

	f.close();
}

int main(int argc, char** argv){
	auto square_ic = [](double x) { return (0 <= x) && (x <= 1) ? 1. : 0.; };
	auto sine_ic = [](double x) { return 1.5 + sin(x); };

	Mesh<AdvectionElement> mesh(50, 0, 2*M_PI, sine_ic);

	create_animation_data(mesh, 0.1, 5, "data.dat");

	return 0;
}