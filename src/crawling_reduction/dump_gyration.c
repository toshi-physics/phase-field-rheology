// dump_gyration.c
// Dump the gyration tensor components of the cells

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "dump.h"
#include "phase_field_model.h"
#include "cell.h"

typedef struct GyrationDump {
  Dump super; // Base struct must be the first element
  bool overwrite;
} GyrationDump;

void gyrationOutput(GyrationDump* dump, PhaseFieldModel* model, int step) {
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
  int clx, cly, count;
  double dx, dy, gxx, gyy, gxy;
  double** cellField;

  fprintf(f, "Cells: %d\n", model->numOfCells);
  fprintf(f, "Timestep: %d\n", step);
  for (int m = 0; m < model->numOfCells; m++) {
    // Assume centre of mass is updated
    cell = model->cells[m];
    clx = cell->lx;
    cly = cell->ly;
    cellField = cell->field[cell->getIndex];
    gxx = 0.0;
    gxy = 0.0;
    gyy = 0.0;
    count = 0;
    for (int i = 0; i < clx; i++) {
      for (int j = 0; j < cly; j++) {
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
    fprintf(f, "%.5f %.5f %.5f\n", gxx, gyy, gxy);
  }
  fclose(f);
  if (dump->overwrite) {
    rename(tmpfile, dump->super.filename);
  }
}

void deleteGyrationDump(GyrationDump* dump) {
  free(dump);
}

DumpFuncs gyrationDumpFuncs =
  {
   .output = (void (*)(Dump*, PhaseFieldModel*, int)) &gyrationOutput,
   .destroy = (void (*)(Dump*)) &deleteGyrationDump
  };

Dump* createGyrationDump(char* filename, int printInc, bool overwrite) {
  GyrationDump* dump = malloc(sizeof *dump);
  setDump(&dump->super, dump, filename, printInc, &gyrationDumpFuncs);
  dump->overwrite = overwrite;
  return (Dump*) dump;
}
