#ifndef MESH_H
#define MESH_H

#include "advection_element.h"
#include "mvector.h"
#include "mmatrix.h"

/*
MeshData holds interpolated values of the element data.
Used for accessing data of the elements outside of the mesh.
*/
struct MeshData {
    MeshData(int n = 0){
        x = vector<double>(n);
        u = vector<double>(n);
    }
    vector<double> x;
    vector<double> u;
};

/*
The Mesh class handles the connectivity of the elements.
Timestepping is done on the mesh, as well as initialisation
of variables for the numerical method.
By templating the class, we can have a mesh for any kind of
element that we like. 
In our case, this gives us a single mesh for AdvectionElement
and BurgerElement.
*/
template<class Element>
class Mesh {

public:
	// Number of elements in the mesh
	int n_elements;
	
    // Empty constructor
    Mesh() : Mesh(0, 0, 0, nullptr){}

    // Constructor initialises the mesh for the given arguments
    Mesh(int mesh_size, double x0, double x1, double (*initial_condition)(double x)){
		n_elements = mesh_size;
        if (mesh_size > 0)
		    create_mesh(mesh_size, x0, x1, initial_condition);
	}
    
	/*
	Calculate interpolated values of the data on the elements.
	Parameter s is a real number between -1 and 1.
	*/
    MeshData get_mesh_data(double s){

        MeshData data(n_elements);
        
		// Loop over all the elements, calculating the interplated values
		// on each element.
        for (int i = 0; i < n_elements; ++i, elements=elements->Right_neighbour_pt){
            data.x[i] = elements->interpolated_x(s);
            data.u[i] = elements->interpolated_u(s);
        }

        return data;
    }

	/*
	Timestep the solution on the mesh.
	To perform a single timestep, there are two loops:

	(1) For each element e
			Calculate the solution on e
	(2) For each element e
			Update the solution on e
	*/
	void timestep(double dt) {

		// Create the inverse of the mass matrix.
		// The matrix is the same for each element since we have
		// a uniform mesh.
		MMatrix Minv(2, 2, 0.);
		Minv(0, 0) = 2; Minv(0, 1) = -1;
		Minv(1, 0) = -1; Minv(1, 1) = 2;
		Minv = Minv * (2. / (elements->X[1] - elements->X[0]));

		// (1) Loop over each element and compute the update values 
		for (int i = 0; i < n_elements; ++i, elements=elements->Right_neighbour_pt){
			
			// Create temporary storage for the updated values
			elements->U_update[0] = elements->U[0];
			elements->U_update[1] = elements->U[1];

			// Flux vector
			MVector F(2);
			F(0) = -0.5; 
			F(1) = 0.5;
			F = F * elements->integrate_flux();

			// Integral of the flux vector
			MVector H(2);
			H(0) = elements->h(elements->Left_neighbour_pt->U[1], elements->U[0]);
			H(1) = -elements->h(elements->U[1], elements->Right_neighbour_pt->U[0]);
			
			// Perform the update step for the temporary storage
			elements->U_update = elements->U_update + dt*Minv*(F+H);
		}

		// (2) Loop over each element and update the values into permanent U
		for (int i = 0; i < n_elements; ++i, elements=elements->Right_neighbour_pt){
			elements->U[0] = elements->U_update[0];
			elements->U[1] = elements->U_update[1];
		}
	}

private:
    // Pointer to the element linked list
    AdvectionElement* elements;

    /*
    Function create_mesh() sets up the connectivity of the mesh elements.
    Since each element is connected by left and right neighbours, as well
    as the start and ended elements being connected, this data structure
    is a circular doubly linked list.
    Elements are added before the "first" element by redirecting the
    neighbour pointers appropriately.
    */
	void create_mesh(int N, double x0, double x1, double (*initial_condition)(double x)) {

		double h = (x1 - x0) / N;

		// Set the element data for the first element
		elements = new Element;
		elements->X[0] = x0;
		elements->X[1] = x0 + h;
		elements->U[0] = initial_condition(elements->X[0]);
		elements->U[1] = initial_condition(elements->X[1]);
		elements->Left_neighbour_pt = elements->Right_neighbour_pt = elements;

		for (int i = 1; i < N; ++i) {
			// Keep track of the last element in the list
			AdvectionElement* last = elements->Left_neighbour_pt;

			// Create the new advection element
			AdvectionElement* e = new Element;
			e->X[0] = x0 + i*h;
			e->X[1] = x0 + (i + 1)*h;
			e->U[0] = initial_condition(e->X[0]);
			e->U[1] = initial_condition(e->X[1]);

			// Update the pointers to maintain the correct sequence of elements
			// - e is the new last element, so left neighbour is last and right neighbour is start
			// - start's new left neighbour is e
			// - last's new right neighbour is e
			e->Right_neighbour_pt = elements;
			e->Left_neighbour_pt = last;
			elements->Left_neighbour_pt = e;
			last->Right_neighbour_pt = e;
		}
	}
};


#endif