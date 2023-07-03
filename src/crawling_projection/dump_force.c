// dump_force.c
// Dump the forces acting on the cells

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "dump.h"
#include "model.h"
#include "cell.h"

typedef struct ForceDump {
  Dump super; // Base struct must be the first element
  bool overwrite;
} ForceDump;

void forceOutput(ForceDump* dump, Model* model, int step) {
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
  int nx, ny;
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
      fprintf(f, "%g %g %d %d ", wx, wy, ix, iy);
    }
    nx = i*2;
    ny = nx+1;
    // Capillary force
    fprintf(f, "%g %g ", model->capillForce[nx], model->capillForce[ny]); 
    // Polar force
    fprintf(f, "%g %g ", model->polarForce[nx], model->polarForce[ny]);
    // Active shear force
    fprintf(f, "%g %g ", model->activeShearForce[nx],
	    model->activeShearForce[ny]);
    // Viscous force
    fprintf(f, "%g %g ", model->viscousForce[nx], model->viscousForce[ny]);
    // Damping force
    fprintf(f, "%g %g\n", model->dampingForce[nx], model->dampingForce[ny]); 
  }
  fclose(f);
  if (dump->overwrite) {
    rename(tmpfile, dump->super.filename);
  }
}

void deleteForceDump(ForceDump* dump) {
  free(dump);
}

DumpFuncs forceDumpFuncs =
  {
   .output = (void (*)(Dump*, Model*, int)) &forceOutput,
   .destroy = (void (*)(Dump*)) &deleteForceDump
  };

Dump* createForceDump(char* filename, int printInc, bool overwrite) {
  ForceDump* dump = malloc(sizeof *dump);
  setDump(&dump->super, dump, filename, printInc, &forceDumpFuncs);
  dump->overwrite = overwrite;
  return (Dump*) dump;
}
