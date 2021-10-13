// dump_deform_field.c
// Dump the deformation tensor field of the cells

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "dump.h"
#include "phase_field_model.h"
#include "cell.h"
#include "array.h"
#include "util.h"

typedef struct DeformFieldDump {
  Dump super; // Base struct must be the first element
  bool overwrite;
} DeformFieldDump;

void deformFieldOutput(DeformFieldDump* dump, PhaseFieldModel* model, 
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
  int clx, cly, cx, cy, x, y;
  int iuu, iu, id, idd, juu, ju, jd, jdd; // Nearest neighbours
  double gradx, grady, gxx, gyy, gxy;
  double** cellField;
  double*** field = create3DDoubleArray(model->lx, model->ly, 3);
  for (int m = 0; m < model->numOfCells; m++) {
    cell = model->cells[m];
    clx = cell->lx;
    cly = cell->ly;
    cx = cell->x;
    cy = cell->y;
    cellField = cell->field[cell->getIndex];
    gxx = 0.0;
    gyy = 0.0;
    gxy = 0.0;
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
	gradx = grad4(i, j, iuu, iu, id, idd, 0, cellField);
	grady = grad4(i, j, juu, ju, jd, jdd, 1, cellField);
	gxx += gradx*gradx;
	gyy += grady*grady;
	gxy += gradx*grady;
      }
    }
    // Add results to total field
    for (int i = 0; i < clx; i++) {
      for (int j = 0; j < cly; j++) {
	x = iwrap(clx, cx+i);
	y = iwrap(cly, cy+j);
	field[x][y][0] += cellField[i][j]*gxx;
	field[x][y][1] += cellField[i][j]*gyy;
	field[x][y][2] += cellField[i][j]*gxy;
      }
    }
  }
  for (int i = 0; i < model->lx; i++) {
    for (int j = 0; j < model->ly; j++) {
      fprintf(f, "%d %d %.5f %.5f %.5f\n", i, j,
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
   .output = (void (*)(Dump*, PhaseFieldModel*, int)) &deformFieldOutput,
   .destroy = (void (*)(Dump*)) &deleteDeformFieldDump
  };

Dump* createDeformFieldDump(char* filename, int printInc, bool overwrite) {
  DeformFieldDump* dump = malloc(sizeof *dump);
  setDump(&dump->super, dump, filename, printInc, &deformFieldDumpFuncs);
  dump->overwrite = overwrite;
  return (Dump*) dump;
}
