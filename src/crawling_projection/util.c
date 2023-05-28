// util.c

#include <stdlib.h>
#include <math.h>
#include "util.h"

inline int iwrap(int len, int i) {
  int remainder = i % len;
  return remainder >= 0 ? remainder : len + remainder;
}

inline int iup(int len, int i) {
  return (i+1 >= len) ? 0 : i+1;
}

inline int idown(int len, int i) {
  return (i-1 < 0) ? len-1 : i-1;
}

inline int idiff(int len, int i1, int i2) {
  int di1 = i1-i2;
  int adi1 = abs(di1);
  int adi2 = len-adi1;
  return adi1 < adi2 ? di1 : -sgni(di1)*adi2;
}

inline double ddiff(double len, double d1, double d2) {
  double dd1 = d1-d2;
  double add1 = fabs(dd1);
  double add2 = len-add1;
  return add1 < add2 ? dd1 : -sgnd(dd1)*add2;
}

inline double laplacian(int i, int j, int iu, int id, int ju, int jd, 
			double** field) {
  return (4.0 * (field[iu][j] + field[id][j] + field[i][ju] + field[i][jd]) +
	  field[iu][ju] + field[iu][jd] + field[id][ju] + field[id][jd] -
	  20.0 * field[i][j]) / 6.0;
}

inline double grad2(int i, int j, int u, int d, int comp, double** field) {
  switch (comp) {
  case 0: return 0.5*(field[u][j] - field[d][j]);
  case 1: return 0.5*(field[i][u] - field[i][d]);
  default: return 0.0;
  }
}

inline double grad4(int i, int j, int uu, int u, int d, int dd, int comp, 
		    double** field) {
  switch (comp) {
  case 0: return (-field[uu][j] + 8.0 * (field[u][j] - field[d][j]) + 
		  field[dd][j]) / 12.0;
  case 1: return (-field[i][uu] + 8.0 * (field[i][u] - field[i][d]) +
		  field[i][dd]) / 12.0;
  default: return 0.0;
  }
}

inline double cgrad4(int i, int j, int uu, int u, int d, int dd, int ic, 
		     int oc, double*** field) {
  switch (oc) {
  case 0: return (-field[uu][j][ic] + 
		  8.0 * (field[u][j][ic] - field[d][j][ic]) + 
		  field[dd][j][ic]) / 12.0;
  case 1: return (-field[i][uu][ic] + 
		  8.0 * (field[i][u][ic] - field[i][d][ic]) +
		  field[i][dd][ic]) / 12.0;
  default: return 0.0;
  }
}

inline double upwind(int i, int j, int uu, int u, int d, int dd, int comp, 
		     double v, double** field) {
  switch (comp) {
  case 0: return v > 0.0 ?
      v * (2.0 * field[u][j] + 3.0 * field[i][j] -
	   6.0 * field[d][j] + field[dd][j]) / 6.0 :
    v * (-field[uu][j] + 6.0 * field[u][j] -
	 3.0 * field[i][j] - 2.0 * field[d][j]) / 6.0;
  case 1: return v > 0.0 ?
      v * (2.0 * field[i][u] + 3.0 * field[i][j] -
	   6.0 * field[i][d] + field[i][dd]) / 6.0 :
    v * (-field[i][uu] + 6.0 * field[i][u] -
	 3.0 * field[i][j] - 2.0 * field[i][d]) / 6.0;
  default: return 0.0;
  }
}

inline int sgni(int val) {
  return (0 < val) - (val < 0);
}

inline double sgnd(double val) {
  return (0.0 < val) - (val < 0.0);
}
