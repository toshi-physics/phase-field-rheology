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

typedef struct VelocityDump {
  Dump super; // Base struct must be the first element
  bool overwrite;
} VelocityDump;

void velocityOutput(VelocityDump* dump, PhaseFieldModel* model, int step) {
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
  int ix, iy;
  double x, y, cx, cy, wx, wy;
  
  fprintf(f, "Cells: %d\n", model->numOfCells);
  fprintf(f, "Timestep: %d\n", step);
  for (int i = 0; i < model->numOfCells; i++) {
    cell = model->cells[i];
    if (dump->overwrite) {
      // Cell CM
      cx = cell->x;
      cy = cell->y;
      x = cx+cell->xcm;
      y = cy+cell->ycm;
      ix = (int) floor(x / model->lx);
      iy = (int) floor(y / model->ly);
      wx = x - ix * model->lx;
      wy = y - iy * model->ly;
      // Output periodic CM and boundary count
      fprintf(f, "%.5f %.5f %d %d ", wx, wy, ix, iy);
    }
    // Advection velocity
    fprintf(f, "%.5e %.5e ", cell->vx, cell->vy);
    // CM velocity
    fprintf(f, "%.5e %.5e\n", cell->drx/model->dt, cell->dry/model->dt);
  }
  fclose(f);
  if (dump->overwrite) {
    rename(tmpfile, dump->super.filename);
  }
}

void deleteVelocityDump(VelocityDump* dump) {
  free(dump);
}

DumpFuncs velocityDumpFuncs =
  {
   .output = (void (*)(Dump*, PhaseFieldModel*, int)) &velocityOutput,
   .destroy = (void (*)(Dump*)) &deleteVelocityDump
  };

Dump* createVelocityDump(char* filename, int printInc, bool overwrite) {
  VelocityDump* dump = malloc(sizeof *dump);
  setDump(&dump->super, dump, filename, printInc, &velocityDumpFuncs);
  dump->overwrite = overwrite;
  return (Dump*) dump;
}
