// linear_lapack.c
// Wrapper of lapack routines

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "linear.h"

// General matrix-vector multiplication
extern void dgemv_(char* trans, int* m, int* n, double* alpha, double* a,
		  int* lda, double* x, int* incx, double* beta, double* y,
		  int* incy);

// General matrix-matrix multiplication
extern void dgemm_(char* transA, char* transB, int* m, int* n, int* k, 
		   double* alpha, double* a, int* lda, double* b, int* ldb,
		   double* beta, double* c, int* ldc);

// General matrix solver
extern void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv,
		  double* b, int* ldb, int* info);

// Helper functions
int* createIntMat(int m, int n);
double* createDoubleMat(int m, int n);

void matvec(double* mat, double* x, double* y, int nrow, int ncol) {
  char trans = 'N';
  double a = 1.0;
  double b = 0.0;
  int one = 1;
  dgemv_(&trans, &nrow, &ncol, &a, mat, &nrow, x, &one, &b, y, &one);
}

void matmat(double* mat1, double* mat2, double* mat3, 
	    int nrow1, int ncol1, int ncol2) {
  ammpbm(mat1, mat2, mat3, 1.0, 0.0, nrow1, ncol1, ncol2);
}

void ammpbm(double* mat1, double* mat2, double* mat3, double a, double b,
	    int nrow1, int ncol1, int ncol2) {
  char trans = 'N';
  dgemm_(&trans, &trans, &nrow1, &ncol2, &ncol1, &a, 
	 mat1, &nrow1, mat2, &ncol1, &b, mat3, &nrow1);
}

void solver(double* mat, double* y, double* x, int nrow, int ncol, int copy) {
  // Work with a copy of the matrix to avoid modifying the original matrix
  double* cmat;
  if (copy) {
    cmat = createDoubleMat(nrow, nrow);
    memcpy(cmat, mat, sizeof(double)*nrow*nrow);
  } else {
    cmat = mat;
  }
  memcpy(x, y, sizeof(double)*nrow*ncol);
  int* ipiv = createIntMat(nrow, 1);
  int info = 0;
  dgesv_(&nrow, &ncol, cmat, &nrow, ipiv, x, &nrow, &info);
  if (copy) {
    free(cmat);
  }
  free(ipiv);
}

void printMatrix(double* a, int m, int n) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      printf("%.5f ", a[m*j+i]);
    }
    printf("\n");
  }
}

inline int* createIntMat(int m, int n) {
  return (int*) malloc(sizeof(int)*m*n);
}

inline double* createDoubleMat(int m, int n) {
  return (double*) malloc(sizeof(double)*m*n);
}
