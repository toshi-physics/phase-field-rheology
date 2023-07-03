// util.h
// Utility math and periodic boundary functions

#ifndef UTIL_H
#define UTIL_H

// For periodic boundaries
int iwrap(int len, int i);
int iup(int len, int i);
int idown(int len, int i);
int idiff(int len, int i1, int i2);
double ddiff(double len, double d1, double d2);

// For taking derivatives
double laplacian(int i, int j, int iu, int id, int ju, int jd, double** field);
double grad2(int i, int j, int u, int d, int comp, double** field);
double grad4(int i, int j, int uu, int u, int d, int dd, int comp, 
	     double** field);
double upwind(int i, int j, int uu, int u, int d, int dd, int comp, 
	      double v, double** field);

int sgni(int val);
double sgnd(double val);

#endif
