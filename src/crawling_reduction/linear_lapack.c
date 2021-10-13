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

void matvec(double* mat, double* x, double* y, int nr) {
  char trans = 'N';
  double alpha = 1.0;
  double beta = 0.0;
  int one = 1;
  dgemv_(&trans, &nr, &nr, &alpha, mat, &nr, x, &one, &beta, y, &one);
}

void matmat(double* mat1, double* mat2, double* mat3, int nr, int nc) {
  char trans = 'N';
  double alpha = 1.0;
  double beta = 0.0;
  dgemm_(&trans, &trans, &nr, &nc, &nr, &alpha,
	 mat1, &nr, mat2, &nr, &beta, mat3, &nr);
}

void ammpbm(double* mat1, double* mat2, double* mat3, double a, double b,
	    int nr, int nc) {
  char trans = 'N';
  dgemm_(&trans, &trans, &nr, &nc, &nr, &a, 
	 mat1, &nr, mat2, &nr, &b, mat3, &nr);
}

void solver(double* mat, double* y, double* x, int nr, int nc, int copy) {
  // Work with a copy of the matrix to avoid modifying the original matrix
  double* cmat;
  if (copy) {
    cmat = createDoubleMat(nr, nr);
    memcpy(cmat, mat, sizeof(double)*nr*nr);
  } else {
    cmat = mat;
  }
  memcpy(x, y, sizeof(double)*nr*nc);
  int* ipiv = createIntMat(nr, 1);
  int info = 0;
  dgesv_(&nr, &nc, cmat, &nr, ipiv, x, &nr, &info);
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
