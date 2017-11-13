#ifndef TESTING_H
#define TESTING_H

#include <iostream>
#include <fstream>
#include "mvector.h"
#include "mmatrix.h"
#include "mbandedmatrix.h"
#include "functions.h"
#include "timing.h"

// Task 3.3.4
void create_timings_laplace_operator_2d(){

    // Output file
    ofstream f("task3.3.4.dat");
    if (!f.is_open()){
        cerr << "Unable to open \"task3.3.4.dat\" for writing.\n"
             << "Exiting...\n";
        exit(-1);
    }

    // Variables for computing cg
    MMatrix A1;
    MBandedMatrix A2;
    MVector b, x0;
    int maxiter = 1000;
    double tol = 1e-6;

    // Timing variables
    double start1 = 0, end1 = 0;
    double start2 = 0, end2 = 0;

    for (int k = 0; k < 40; ++k ){
        int n = 5 * (k+1);
        A1 = laplacian_2d(n);
        A2 = laplacian_2d_banded(n);
        b = create_vector1(n * n);
        x0 = MVector(n * n, 0.);

        // Time CG using MMatrix
        start1 = Timer();
        VectorSolution sol1 = cg(A1, b, x0, maxiter, tol);
        end1 = Timer();

        // Time CG using MBandedMatrix
        start2 = Timer();
        VectorSolution sol2 = cg(A2, b, x0, maxiter, tol);
        end2 = Timer();

        // Output format: tab separated values
        //  size of matrix;
        //  time to execute cg for MMatrix;
        //  iterations to converge for MMatrix;
        //  time to execute cg for MBandedMatrix;
        //  iterations to converge for MBandedMatrix;
        f << n << "\t" << (end1 - start1) << "\t" << sol1.iteration_count << "\t"
                       << (end2 - start2) << "\t" << sol2.iteration_count << endl;
    }

    f.close();
}

// Task 3.3.5
void create_plot_data_laplace_operator_2d(){

    // Output file
    ofstream f("task3.3.5.dat");
    if (!f.is_open()){
        cerr << "Unable to open \"task3.3.5.dat\" for writing.\n"
             << "Exiting...\n";
        exit(-1);
    }

    // Variables for computing cg
    int n = 25;
    MBandedMatrix A = laplacian_2d_banded(n);
    MVector b = create_vector1(n * n);
    MVector x0(n * n, 0.);
    int maxiter = 1000;
    double tol = 1e-6;

    VectorSolution sol = cg(A, b, x0, maxiter, tol);

    for (int i = 0; i < n; ++i){
        for (int j = 0; j < n - 1; ++j){
            f << sol.solution(i * n + j) << "\t";
        }
        f << sol.solution(i * n + n - 1) << endl;
    }

    f.close();
}

#endif