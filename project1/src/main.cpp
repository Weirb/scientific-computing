#include <iostream>
#include "mvector.h"
#include "mmatrix.h"
#include "mbandedmatrix.h"
#include "functions.h"
#include "testing.h"

using namespace std;

void generate_all_data() {
	// This is memory intensive, may cause the machine to hold or crash.
	// Visual Studio can handle this with 8GB of RAM
	
	try {
		cout << "create_cg_plot_laplacian_1d" << endl;
		create_cg_plot_laplacian_1d();

		cout << "create_cg_table_laplacian_1d" << endl;
		create_cg_table_laplacian_1d();

		cout << "create_table_laplacian_1d_banded" << endl; 
		create_table_laplacian_1d_banded();

		cout << "create_timings_laplace_operator_2d" << endl;
		create_timings_laplace_operator_2d();

		cout << "create_plot_data_laplace_operator_2d" << endl; 
		create_plot_data_laplace_operator_2d();

		cout << "create_cg_plot_matrix2" << endl;
		create_cg_plot_matrix2();

		cout << "create_cg_table_matrix2" << endl;
		create_cg_table_matrix2();

	} catch(const bad_alloc& e){
		cout << e.what() 
			 << "Memory error, could not generate data. Exiting..." << endl;
	}
}

int main(int argc, char** argv){

	generate_all_data();

	return 0;
}