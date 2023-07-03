// linear.h
// A collection of methods for doing linear algebra
// Matrix A should be in column major format and in squared form

#ifndef LINEAR_H
#define LINEAR_H

// Compute y := mat*x
void matvec(double* mat, double* x, double* y, int nrow, int ncol);

// Compute mat3 := mat1*mat2
void matmat(double* mat1, double* mat2, double* mat3,
	    int nrow1, int ncol1, int ncol2);

// Compute mat3 := a*mat1*mat2+b*mat3
void ammpbm(double* mat1, double* mat2, double* mat3, double a, double b,
	    int nrow1, int ncol1, int ncol2);

// Solve the system y = mat*x for x
void solver(double* mat, double* y, double* x, int nrow, int ncol, int copy);

// Print a matrix of dimensions m x n (m = number of rows, n = number of cols)
void printMatrix(double* a, int nrow, int ncol);

#endif
