#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "mvector.h"
#include "mmatrix.h"
#include "mbandedmatrix.h"


MMatrix laplacian_1d(int n);
MBandedMatrix laplacian_1d_banded(int n);

MMatrix create_matrix2(int n, double m);

MMatrix laplacian_2d(int n);
MBandedMatrix laplacian_2d_banded(int n);

MVector create_vector1(int n);
MVector create_vector2(int n);

MMatrix kronecker(MMatrix A, MMatrix B);

#endif
