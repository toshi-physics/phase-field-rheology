// dump_velocity_field.c
// Dump the velocity field of the cells

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

typedef struct VelocityFieldDump {
  Dump super; // Base struct must be the first element
  bool overwrite;
} VelocityFieldDump;

void velocityFieldOutput(VelocityFieldDump* dump, Model* model, 
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
  // Sum over all the phase fields, each weighted by the velocity vector
  // of the cell  
  Cell* cell;
  int clx, cly, cx, cy, x, y;
  double cvx, cvy, cvcmx, cvcmy;
  double** cellField;
  double*** field = create3DDoubleArray(model->lx, model->ly, 5);
  for (int m = 0; m < model->numOfCells; m++) {
    cell = model->cells[m];
    clx = cell->lx;
    cly = cell->ly;
    cx = cell->x;
    cy = cell->y;
    cvx = cell->vx;
    cvy = cell->vy;
    cvcmx = cell->drx/model->dt;
    cvcmy = cell->dry/model->dt;
    cellField = cell->field[cell->getIndex];
    for  (int i = 0; i < clx; i++) {
      for (int j = 0; j < cly; j++) {
	x = iwrap(model->lx, cx+i);
	y = iwrap(model->ly, cy+j);
	field[x][y][0] += cellField[i][j];
	field[x][y][1] += cellField[i][j] * cvx;
	field[x][y][2] += cellField[i][j] * cvy;
	field[x][y][3] += cellField[i][j] * cvcmx;
	field[x][y][4] += cellField[i][j] * cvcmy;
	
      }
    }
  }
  for (int i = 0; i < model->lx; i++) {
    for (int j = 0; j < model->ly; j++) {
      fprintf(f, "%d %d %g %g %g %g\n", i, j, 
	      field[i][j][1] / field[i][j][0],
	      field[i][j][2] / field[i][j][0],
	      field[i][j][3] / field[i][j][0],
	      field[i][j][4] / field[i][j][0]);
    }
    fprintf(f,"\n");
  }
  fclose(f);
  free(field);
  if (dump->overwrite) {
    rename(tmpfile, dump->super.filename);
  }
}

void deleteVelocityFieldDump(VelocityFieldDump* dump) {
  free(dump);
}

DumpFuncs velocityFieldDumpFuncs =
  {
   .output = (void (*)(Dump*, Model*, int)) &velocityFieldOutput,
   .destroy = (void (*)(Dump*)) &deleteVelocityFieldDump
  };

Dump* createVelocityFieldDump(char* filename, int printInc, bool overwrite) {
  VelocityFieldDump* dump = malloc(sizeof *dump);
  setDump(&dump->super, dump, filename, printInc, &velocityFieldDumpFuncs);
  dump->overwrite = overwrite;
  return (Dump*) dump;
}
