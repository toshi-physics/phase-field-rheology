// dump_deform.c
// Dump the deformation tensor of each cell

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "dump.h"
#include "model.h"
#include "cell.h"
#include "util.h"

typedef struct DeformDump {
  Dump super; // Base struct must be the first element
  bool overwrite;
} DeformDump;

void deformOutput(DeformDump* dump, Model* model, int step) {
  char tmpfile [PF_DIR_SIZE];
  FILE* f;
  if (dump->overwrite) {
    strcpy(tmpfile, dump->super.filename);
    strcat(tmpfile, ".tmp");
    f = fopen(tmpfile, "w");
  } else {
    f = fopen(dump->super.filename, "a");
  }
  
  Cell* cell;
  int clx, cly;
  int iuu, iu, id, idd, juu, ju, jd, jdd; // Nearest neighbours
  double gphix, gphiy, sxx, sxy, s, s1, s1x, s1y, s2, s2x, s2y;
  double** cellField;

  fprintf(f, "Cells: %d\n", model->numOfCells);
  fprintf(f, "Timestep: %d\n", step);
  for (int m = 0; m < model->numOfCells; m++) {
    cell = model->cells[m];
    clx = cell->lx;
    cly = cell->ly;
    cellField = cell->field[cell->getIndex];
    sxx = 0.0;
    sxy = 0.0;
    for (int i = 2; i < clx-2; i++) {
      iu = iup(clx, i);
      iuu = iup(clx, iu);
      id = idown(clx, i);
      idd = idown(clx, id);
      for (int j = 2; j < cly-2; j++) {
	ju = iup(cly, j);
	juu = iup(cly, ju);
	jd = idown(cly, j);
	jdd = idown(cly, jd);
	gphix = grad4(i, j, iuu, iu, id, idd, 0, cellField);
	gphiy = grad4(i, j, juu, ju, jd, jdd, 1, cellField);
	sxx += gphiy * gphiy - gphix * gphix;
	sxy -= gphix * gphiy;
      }
    }
    sxx *= 0.5;
    // Compute eigenvalues and eigenvectors
    s = sqrt(sxx * sxx + sxy * sxy);
    s1y = (s - sxx) / sxy;
    s1 = sqrt(1.0 + s1y * s1y); // s1x = 1.0
    s1x = 1.0 / s1;
    s1y /= s1;
    s2y = -(s + sxx) / sxy;
    s2 = sqrt(1.0 + s2y * s2y); // s2x = 1.0
    s2x = 1.0 / s2;
    s2y /= s2;
    fprintf(f, "%f %f %f %f %f %f %f\n", sxx, sxy, s, s1x, s1y, s2x, s2y);
  }
  fclose(f);
  if (dump->overwrite) {
    rename(tmpfile, dump->super.filename);
  }
}

void deleteDeformDump(DeformDump* dump) {
  free(dump);
}

DumpFuncs deformDumpFuncs =
  {
   .output = (void (*)(Dump*, Model*, int)) &deformOutput,
   .destroy = (void (*)(Dump*)) &deleteDeformDump
  };

Dump* createDeformDump(char* filename, int printInc, bool overwrite) {
  DeformDump* dump = malloc(sizeof *dump);
  setDump(&dump->super, dump, filename, printInc, &deformDumpFuncs);
  dump->overwrite = overwrite;
  return (Dump*) dump;
}
