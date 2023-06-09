// dump_gyration_field.c
// Dump the gyration tensor field of the cells

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "dump.h"
#include "model.h"
#include "cell.h"
#include "array.h"
#include "util.h"

typedef struct GyrationFieldDump {
  Dump super; // Base struct must be the first element
  bool overwrite;
} GyrationFieldDump;

void gyrationFieldOutput(GyrationFieldDump* dump, Model* model, 
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
  // Sum over all the phase fields, each weighted by the gyration tensor
  // of the cell  
  Cell* cell;
  int clx, cly, cx, cy, x, y, count;
  double dx, dy, gxx, gyy, gxy;
  double** cellField;
  double*** field = create3DDoubleArray(model->lx, model->ly, 3);
  for (int m = 0; m < model->numOfCells; m++) {
    // Assume centre of mass is updated
    cell = model->cells[m];
    clx = cell->lx;
    cly = cell->ly;
    cx = cell->x;
    cy = cell->y;
    cellField = cell->field[cell->getIndex];
    gxx = 0.0;
    gxy = 0.0;
    gyy = 0.0;
    count = 0;
    for (int i = 2; i < clx-2; i++) {
      for (int j = 2; j < cly-2; j++) {
	if (cellField[i][j] > cell->incell) {
	  dx = i+0.5-cell->xcm;
	  dy = j+0.5-cell->ycm;
	  gxx += dx*dx;
	  gyy += dy*dy;
	  gxy += dx*dy;
	  count++;
	}
      }
    }
    gxx /= (double) count;
    gxy /= (double) count;
    gyy /= (double) count;
    for (int i = 2; i < clx-2; i++) {
      for (int j = 2; j < cly-2; j++) {
	x = iwrap(model->lx, cx+i);
	y = iwrap(model->ly, cy+j);
	field[x][y][0] += cellField[i][j] * gxx;
	field[x][y][1] += cellField[i][j] * gyy;
	field[x][y][2] += cellField[i][j] * gxy;	
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

void deleteGyrationFieldDump(GyrationFieldDump* dump) {
  free(dump);
}

DumpFuncs gyrationFieldDumpFuncs =
  {
   .output = (void (*)(Dump*, Model*, int)) &gyrationFieldOutput,
   .destroy = (void (*)(Dump*)) &deleteGyrationFieldDump
  };

Dump* createGyrationFieldDump(char* filename, int printInc, bool overwrite) {
  GyrationFieldDump* dump = malloc(sizeof *dump);
  setDump(&dump->super, dump, filename, printInc, &gyrationFieldDumpFuncs);
  dump->overwrite = overwrite;
  return (Dump*) dump;
}
