#ifndef TESTING_H
#define TESTING_H

#include <iostream>
#include <fstream>
#include "mvector.h"
#include "mmatrix.h"
#include "mbandedmatrix.h"
#include "functions.h"
#include "timing.h"
#include <vector>
#include <string>

// Task 3.2.1.3
void test_vector_overloads(){
/*
Test the addition and subtraction of MVectors.
*/
	MVector v(3, 0.);
    v(0) = 0.1;
    v(1) = 4.8;
    v(2) = 3.7;

    MVector w(3);
    w(0) = 3.1;
    w(1) = 8.5;
    w(2) = 3.6;

    MVector x(3);
    x(0) = 5.8;
    x(1) = 7.4;
    x(2) = 12.4;

    MVector u = 4.7*v + 1.3*w - 6.7*x;

    double tol = 1e-5;
    bool u1 = u(1) - 4.7 * 0.1 + 1.3 * 3.1 - 6.7 * 5.8 < tol;
    bool u2 = u(2) - 4.7 * 4.8 + 1.3 * 8.5 - 6.7 * 7.4 < tol;
    bool u3 = u(3) - 4.7 * 3.7 + 1.3 * 3.6 - 6.7 * 12.4 < tol;

    if (!(u1 && u2 && u3)){
        cout << "Error in addition and subtraction of vectors.\n"
             << "Implementation is incorrect.\n";
    }
}

// Task 3.2.2.4
void test_norms_dotproduct(){
/*
Test the norms and dot product of MVectors
*/
	MVector u(3);
    u(0) = 1.5;
    u(1) = 1.3;
    u(2) = 2.8;

    MVector v(3);
    v(0) = 6.5;
    v(1) = 2.7;
    v(2) = 2.9;

    MVector w(3);
    w(0) = 0.1;
    w(1) = -7.2;
    w(2) = 3.4;

    // Answers correct to 6sf
    double tol = 1e-5;

    // Test u
    cout << "Test of L2Norm(u): ";
    if (u.L2Norm() - 3.4322 < tol)
        cout << "correct";
    else
        cout << "incorrect";
    cout << endl;

    // Test v
    cout << "Test of L2Norm(v): ";
    if (v.L2Norm() - 7.61249 < tol)
        cout << "correct";
    else
        cout << "incorrect";
    cout << endl;

    // Test w
    cout << "Test of L2Norm(w): ";
    if (w.L2Norm() - 7.96304 < tol)
        cout << "correct";
    else
        cout << "incorrect";
    cout << endl;

    // Test value of alpha
    cout << "Test dot product: ";
    double alpha = dot(u, u) / dot(v, w);
    if (alpha - (-1.31915) < tol)
        cout << "correct";
    else
        cout << "incorrect";
    cout << endl;
}

// Task 3.2.3.3
void test_matrix_vector_product() {
/* MATLAB test code
m = 4; n = 3;
A = zeros(m, n);
x = [0.5, 1.6, 3.2]';
for i = 1:m
    for j = 1:n
        A(i,j) = 3*(i-1) + (j-1);
    end
end
A*x
*/
    int m = 4, n = 3;
    MMatrix A(m, n);

    for (int i = 0; i < m; ++i){
        for (int j = 0; j < n; ++j){
            A(i, j) = 3 * i + j;
        }
    }

    MVector x(3);
    x(0) = 0.5;
    x(1) = 1.6;
    x(2) = 3.2;

    MVector b = A*x;
}

// Task 3.2.4.5
void test_cg_algorithm(){

    int n = 25;
    MMatrix A = laplacian_1d(n);
    MVector b = create_vector1(n);
    MVector x0(n, 0.);
    double tol = 1e-6;

    VectorSolution sol = cg(A, b, x0, 1000, tol);

    if (sol.iteration_count == 12){
        cout << "Conjugate gradient took 12 iterations to converge.\n"
             << "Algorithm converged in expected manner.\n";
    } else {
        cout << "Conjugate gradient algorithm did "
                 "not take 12 iterations to converge.\n";
    }
}

// Task 3.2.4.6
void create_cg_plot_laplacian_1d(){

    // CG variables
    MMatrix A;
    MVector b, x0;
    double tol = 1e-6;
    VectorSolution sol;

    // Output files
    ofstream f("task3.2.4.6.dat");
    if (!f.is_open()){
        cerr << "Unable to open \"task3.2.4.6.dat\" for writing.\n"
             << "Exiting...\n";
        exit(-1);
    }

    // Size of problems to solve
    vector<int> ns = {10, 25, 100};

    for (size_t k = 0; k < ns.size(); ++k){

        int n = ns[k];
        A = laplacian_1d(n);
        b = create_vector1(n);
        x0 = MVector(n, 0.);
        
        // Compute solution
        sol = cg(A, b, x0, 1000, tol);

        f << n << endl;
        for (int i = 0; i < n; ++ i){
            // Output format: tab separated values
            //  i;
            //  (i+1)/(n+1);
            //  x_i;
            f << i << "\t"
              << (i + 1) / double(n + 1) << "\t"
              << sol.solution(i) << endl;
        }
        
    }
    f.close();
}

// Task 3.2.4.7
void create_cg_table_laplacian_1d(){

    // CG variables
    MMatrix A;
    MVector b, x0;
    double tol = 1e-6;
    VectorSolution sol;

    double start, end;

    // Output files
    ofstream f("task3.2.4.7.dat");
    if (!f.is_open()){
        cerr << "Unable to open \"task3.2.4.7.dat\" for writing.\n"
             << "Exiting...\n";
        exit(-1);
    }

    for (int k = 0; k < 15; ++k){

        int n = 5 * (k+1);
        A = laplacian_1d(n);
        b = create_vector1(n);
        x0 = MVector(n, 0.);
        
        // Compute solution
        start = Timer();
        sol = cg(A, b, x0, 1000, tol);
        end = Timer();

        // Output format: tab separated values
        //  n;
        //  iterations to converge;
        //  time to converge;
        f << n << "\t"
          << sol.iteration_count << "\t"
          << (end-start) << endl;
    }
    f.close();
}

// Task 3.2.4.8, part 1
void create_cg_plot_matrix2(){

    // CG variables
    MMatrix A;
    MVector b, x0;
    double tol = 1e-6;
    VectorSolution sol;

    // Output files
    ofstream f("task3.2.4.6.dat");
    if (!f.is_open()){
        cerr << "Unable to open \"task3.2.4.6.dat\" for writing.\n"
             << "Exiting...\n";
        exit(-1);
    }

    // Size of problems to solve
    vector<int> ns = {10, 25, 100};

    for (size_t k = 0; k < ns.size(); ++k){

        int n = ns[k];
        A = laplacian_1d(n);
        b = create_vector1(n);
        x0 = MVector(n, 0.);
        
        // Compute solution
        sol = cg(A, b, x0, 1000, tol);

        f << n << endl;
        for (int i = 0; i < n; ++ i){
            // Output format: tab separated values
            //  i;
            //  (i+1)/(n+1);
            //  x_i;
            f << i << "\t"
              << (i + 1) / double(n + 1) << "\t"
              << sol.solution(i) << endl;
        }
        
    }
    f.close();
}

// Task 3.2.4.8, part 2
void create_cg_table_matrix2(){
    
    int n = 100;

    // CG variables
    MMatrix A;
    MVector b = create_vector2(n);
    MVector x0(n, 0.);
    double tol = 1e-3;
    VectorSolution sol;

    double start, end;
    double m;

    // Output files
    ofstream f("task3.2.4.8-2.dat");
    if (!f.is_open()){
        cerr << "Unable to open \"task3.2.4.8-2.dat\" for writing.\n"
             << "Exiting...\n";
        exit(-1);
    }
    
    for (int k = 0; k < 10; ++k){

        // m = 100 * pow(10,-k);
        // m = 10 - k;
        // m = 5 - 0.1*k;
        // m = 4.8 - 0.01*k;
        // m = 4.75 - 0.001*k;
        // m = 4.743 - 0.0001*k;
        m = 4.7429 - 0.00001*k;
        A = create_matrix2(n, m);
        
        // Compute solution
        start = Timer();
        sol = cg(A, b, x0, 1000000, tol);
        end = Timer();

        // Output format: tab separated values
        //  n;
        //  m;
        //  iterations to converge;
        //  time to converge;
        f << n << "\t"
          << m << "\t"
          << sol.iteration_count << "\t"
          << (end-start) << endl;
    }
    f.close();
}

// Task 3.3.4
void create_timings_laplace_operator_2d(){
// Find execution times for computing the solution using the CG algorithm
// for MMatrix and MBandedMatrix

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
        f << n << "\t" << (end1 - start1) << "\t" 
                       << sol1.iteration_count << "\t"
                       << (end2 - start2) << "\t" 
                       << sol2.iteration_count << endl;
    }

    f.close();
}

// Task 3.3.5
void create_plot_data_laplace_operator_2d(){
// Generate a plot of the solution to the 2d laplace equation whose right 
// hand side is the vector whose elements are 1/(n+1)^2

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