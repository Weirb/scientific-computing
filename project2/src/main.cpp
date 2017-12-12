#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include "mvector.h"
#include "mmatrix.h"
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

	int N = 200;
	double dx = 2 * M_PI / (double)N;
	// double dt = 1e-4;
	double dt = dx / 16;

	double T = 4;
	int iteration_count = T / dt;

	Mesh<AdvectionElement> mesh1(N, 0, 2*M_PI, sine_ic);
	Mesh<BurgerElement> mesh2(N, 0, 2*M_PI, sine_ic);

	// Mesh<AdvectionElement> mesh2(10, 0, 2*M_PI, sine_ic);
	// MeshData data = mesh2.get_mesh_data(0);
	// for (int i = 0; i < 10; ++i){
	// 	cout << data.x[i] << ", " << data.u[i] << endl;
	// }

	create_animation_data(mesh1, dt, iteration_count, "data1.dat");
	create_animation_data(mesh2, dt, iteration_count, "data2.dat");

	return 0;
}