// dump_deform_field.c
// Dump the deformation tensor field of the cells

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "dump.h"
#include "model.h"
#include "cell.h"
#include "array.h"
#include "util.h"

typedef struct DeformFieldDump {
  Dump super; // Base struct must be the first element
  bool overwrite;
} DeformFieldDump;

void deformFieldOutput(DeformFieldDump* dump, Model* model, 
		       int step) {
  char tmpfile [PF_DIR_SIZE];
  if (dump->overwrite) {
    strcpy(tmpfile, dump->super.filename);
    strcat(tmpfile, ".tmp");
  } else {
    strcpy(tmpfile, dump->super.filename);
    char suffix [80];
    sprintf(suffix, ".%d", step);
    strcat(tmpfile, suffix);
  }
  FILE* f = fopen(tmpfile, "w");
  // Sum over all the phase fields, each weighted by the deformation tensor
  // of the cell
  Cell* cell;
  int clx, cly, cx, cy, x, y, buf;
  int iuu, iu, id, idd, juu, ju, jd, jdd; // Nearest neighbours
  double gphix, gphiy, sxx, syy, sxy;
  double** cellField;
  double*** field = create3DDoubleArray(model->lx, model->ly, 3);
  for (int m = 0; m < model->numOfCells; m++) {
    cell = model->cells[m];
    clx = cell->lx;
    cly = cell->ly;
    cx = cell->x;
    cy = cell->y;
    buf = cell->haloWidth;
    cellField = cell->field[cell->getIndex];
    sxx = 0.0;
    syy = 0.0;
    sxy = 0.0;
    for (int i = buf; i < clx-buf; i++) {
      iu = iup(clx, i);
      iuu = iup(clx, iu);
      id = idown(clx, i);
      idd = idown(clx, id);
      for (int j = buf; j < cly-buf; j++) {
	ju = iup(cly, j);
	juu = iup(cly, ju);
	jd = idown(cly, j);
	jdd = idown(cly, jd);
	gphix = grad4(i, j, iuu, iu, id, idd, 0, cellField);
	gphiy = grad4(i, j, juu, ju, jd, jdd, 1, cellField);
	sxx += gphix * gphix; // Use +ve sign here to save storage space
	syy += gphiy * gphiy;
	sxy += gphix * gphiy;
      }
    }    
    // Add results to total field
    for (int i = buf; i < clx-buf; i++) {
      for (int j = buf; j < cly-buf; j++) {
	x = iwrap(model->lx, cx+i);
	y = iwrap(model->ly, cy+j);
	field[x][y][0] += cellField[i][j] * sxx;
	field[x][y][1] += cellField[i][j] * syy;
	field[x][y][2] += cellField[i][j] * sxy;
      }
    }
  }
  for (int i = 0; i < model->lx; i++) {
    for (int j = 0; j < model->ly; j++) {
      fprintf(f, "%d %d %g %g %g\n", i, j,
	      field[i][j][0], field[i][j][1], field[i][j][2]);
    }
    fprintf(f,"\n");
  }
  fclose(f);
  free(field);
  if (dump->overwrite) {
    rename(tmpfile, dump->super.filename);
  }
}

void deleteDeformFieldDump(DeformFieldDump* dump) {
  free(dump);
}

DumpFuncs deformFieldDumpFuncs =
  {
   .output = (void (*)(Dump*, Model*, int)) &deformFieldOutput,
   .destroy = (void (*)(Dump*)) &deleteDeformFieldDump
  };

Dump* createDeformFieldDump(char* filename, int printInc, bool overwrite) {
  DeformFieldDump* dump = malloc(sizeof *dump);
  setDump(&dump->super, dump, filename, printInc, &deformFieldDumpFuncs);
  dump->overwrite = overwrite;
  return (Dump*) dump;
}
