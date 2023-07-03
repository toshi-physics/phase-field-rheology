// dump_velocity.c
// Dump the velocities (reduced forces) of the cells

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "dump.h"
#include "phase_field_model.h"
#include "cell.h"
#include "util.h"

typedef struct DeformDump {
  Dump super; // Base struct must be the first element
  bool overwrite;
} DeformDump;

void deformOutput(DeformDump* dump, PhaseFieldModel* model, int step) {
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
  double gradx, grady, gxx, gyy, gxy;
  double** cellField;

  fprintf(f, "Cells: %d\n", model->numOfCells);
  fprintf(f, "Timestep: %d\n", step);
  for (int m = 0; m < model->numOfCells; m++) {
    cell = model->cells[m];
    clx = cell->lx;
    cly = cell->ly;
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
    fprintf(f, "%.5f %.5f %.5f\n", gxx, gyy, gxy);
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
   .output = (void (*)(Dump*, PhaseFieldModel*, int)) &deformOutput,
   .destroy = (void (*)(Dump*)) &deleteDeformDump
  };

Dump* createDeformDump(char* filename, int printInc, bool overwrite) {
  DeformDump* dump = malloc(sizeof *dump);
  setDump(&dump->super, dump, filename, printInc, &deformDumpFuncs);
  dump->overwrite = overwrite;
  return (Dump*) dump;
}
