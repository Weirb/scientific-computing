#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include "mvector.h"
#include "element.h"
#include "mesh.h"

using namespace std;

/*
Create a file containing the solution at time T=N*dt.
*/
template<class Element>
void create_data_file(Mesh<Element> &mesh, double dt, int N, string filename) {

	ofstream f(filename);
	if (!f.is_open()) {
		cout << "Could not open file: " << filename << endl;
		exit(-1);
	}

	// Timestep the mesh, calculating the solution at the desired time
	for (int i = 0; i < N; ++i)
		mesh.timestep(dt);

	// Write the solution data to the file
	MeshData data = mesh.get_mesh_data(0);
	for (int i = 0; i < mesh.n_elements; ++i) {
		f << data.x[i] << "," << data.u[i]<< endl;
	}

	f.close();
}

/*
Create all of the desired files with the solution data.
We want to solve the advection and Burgers' equation for
multiple initial conditions and at several times.
Perform all required generation tasks here.
*/
void generate_data(){
	// Initial condition lambdas
	auto square_ic = [](double x) { return (0 <= x) && (x <= 1) ? 1. : 0.; };
	auto sine_ic = [](double x) { return 1.5 + sin(x); };
	
	// Solution times and filenames
	vector<double> t1 = {0., 0.25, 0.5, 1.};
	vector<double> t2 = {0., 0.25, 0.5, 1., 1.25, 1.5, 1.75, 2.};
	vector<string> files_a = {	"data/advec_sin_t_0.dat", "data/advec_sin_t_1.dat", 
								"data/advec_sin_t_2.dat", "data/advec_sin_t_3.dat",
							  	"data/advec_sq_t_0.dat", "data/advec_sq_t_1.dat", 
								"data/advec_sq_t_2.dat", "data/advec_sq_t_3.dat"};
	vector<string> files_b = {	"data/burgers_sin_t_0.dat", "data/burgers_sin_t_1.dat", 
								"data/burgers_sin_t_2.dat", "data/burgers_sin_t_3.dat",
								"data/burgers_sin_t_4.dat", "data/burgers_sin_t_5.dat", 
								"data/burgers_sin_t_6.dat", "data/burgers_sin_t_7.dat",
								"data/burgers_sq_t_0.dat", "data/burgers_sq_t_1.dat", 
								"data/burgers_sq_t_2.dat", "data/burgers_sq_t_3.dat",
								"data/burgers_sq_t_4.dat", "data/burgers_sq_t_5.dat", 
								"data/burgers_sq_t_6.dat", "data/burgers_sq_t_7.dat"};
	
	// Solution properties
	int N = 100;
	double dx = 2 * M_PI / (double)N;
	double dt = 1e-4;
	cout << (dt / dx) << endl;

	Mesh<AdvectionElement> mesh_a;
	Mesh<BurgersElement> mesh_b;

	// Create first set for advection element data
	for (int i = 0; i < t1.size(); ++ i) {
		int iteration_count = t1[i] / dt;
		mesh_a = Mesh<AdvectionElement>(N, 0, 2 * M_PI, sine_ic);
		create_data_file(mesh_a, dt, iteration_count, files_a[i]);
	}
	
	for (int i = 0; i < t1.size(); ++ i) {
		int iteration_count = t1[i] / dt;
		mesh_a = Mesh<AdvectionElement>(N, 0, 2 * M_PI, square_ic);
		create_data_file(mesh_a, dt, iteration_count, files_a[t1.size()+i]);
	}

	// Create second set for Burgers' element data
	for (int i = 0; i < t2.size(); ++ i) {
		int iteration_count = t2[i] / dt;
		mesh_b = Mesh<BurgersElement>(N, 0, 2 * M_PI, sine_ic);
		create_data_file(mesh_b, dt, iteration_count, files_b[i]);
	}

	for (int i = 0; i < t2.size(); ++ i) {
		int iteration_count = t2[i] / dt;
		mesh_b = Mesh<BurgersElement>(N, 0, 2 * M_PI, square_ic);
		create_data_file(mesh_b, dt, iteration_count, files_b[t2.size()+i]);
	}
}

int main(int argc, char** argv){

	generate_data();
	return 0;
}