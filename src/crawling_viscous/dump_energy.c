// dump_energy.c
// Dump the cm of the cells

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "array.h"
#include "dump.h"
#include "model.h"
#include "cell.h"
#include "util.h"
#include "constant.h"

typedef struct EnergyDump {
  Dump super; // Base struct must be the first element
  int lx;
  int ly;
  double** phi2Field; // Store the sum of phi^2 field
  double*** gphiField; // Store the sum of grad phix and grad phiy
} EnergyDump;

void energyOutput(EnergyDump* dump, Model* model, int step) {
  FILE* f;
  f = fopen(dump->super.filename, "a");
  
  // Reset fields
  for (int i = 0; i < dump->lx; i++) {
    for (int j = 0; j < dump->ly; j++) {
      dump->phi2Field[i][j] = 0.0;
    }
  }
  
  // Compute the auxillary field (sum of phi^2)
  Cell* cell;
  double phi, phi2, gphix, gphiy;
  int clx, cly, x, y, cx, cy;
  int iu, iuu, id, idd, ju, juu, jd, jdd;
  double** cellField;
  for (int n = 0; n < model->numOfCells; n++) {
    cell = model->cells[n];
    clx = cell->lx;
    cly = cell->ly;
    cx = cell->x;
    cy = cell->y;
    cellField = cell->field[cell->getIndex];
    for (int i = 0; i < clx; i++) {
      iu = iup(clx, i);
      iuu = iup(clx, iu);
      id = idown(clx, i);
      idd = idown(clx, id);
      x = iwrap(model->lx, cx+i);
      for (int j = 0; j < cly; j++) {
	ju = iup(cly, j);
	juu = iup(cly, ju);
	jd = idown(cly, j);
	jdd = idown(cly, jd);
	y = iwrap(model->ly, cy+j);
	phi = cellField[i][j];
	phi2 = phi*phi;
	gphix = grad4(i, j, iuu, iu, id, idd, 0, cellField);
	gphiy = grad4(i, j, juu, ju, jd, jdd, 1, cellField);
	dump->phi2Field[x][y] += phi2;
	dump->gphiField[x][y][0] += gphix;
	dump->gphiField[x][y][1] += gphiy;
      } // Close loop over j
    } // Close loop over i
  } // Close loop over cell m

  // Compute the free energy
  double dphi, dphi2, gphi2, cellVolume, piR2, dV;
  double cahnHilliard = 0.0;
  double volumeConstraint = 0.0;
  double repulsion = 0.0;
  double adhesion = 0.0;
  double thickness2 = model->thickness * model->thickness;
  double phi4Sum = 0.0;
  double gphi2Sum = 0.0;
  for (int n = 0; n < model->numOfCells; n++) {
    cell = model->cells[n];
    clx = cell->lx;
    cly = cell->ly;
    cx = cell->x;
    cy = cell->y;
    cellField = cell->field[cell->getIndex];
    cellVolume = 0.0;
    piR2 = PF_PI * cell->radius * cell->radius;
    for (int i = 1; i < clx-1; i++) {
      iu = iup(clx, i);
      iuu = iup(clx, iu);
      id = idown(clx, i);
      idd = idown(clx, id);
      x = iwrap(model->lx, cx+i);
      for (int j = 1; j < cly-1; j++) {
	ju = iup(cly, j);
	juu = iup(cly, ju);
	jd = idown(cly, j);
	jdd = idown(cly, jd);
	y = iwrap(model->ly, cy+j);
	
	phi = cellField[i][j];
	phi2 = phi * phi;
	dphi = phi - 1.0;
	dphi2 = dphi * dphi;
	phi4Sum += phi2 * phi2;
	gphix = grad4(i, j, iuu, iu, id, idd, 0, cellField);
	gphiy = grad4(i, j, juu, ju, jd, jdd, 1, cellField);
	gphi2 = gphix * gphix + gphiy * gphiy;
	gphi2Sum += gphi2;

	// Cahn-Hilliard term
	cahnHilliard += phi2 * dphi2 + thickness2 * gphi2;
	
	// Calculate cell volume
	cellVolume += phi2;
	
	// Repulsion term
	repulsion += phi2 * dump->phi2Field[x][y];

	// Adhesion term
	adhesion += gphix * dump->gphiField[x][y][0] + 
	  gphiy * dump->gphiField[x][y][1];
      } // Close loop over j
    } // Close loop over i
    dV = 1.0 - cellVolume / piR2;
    volumeConstraint += piR2 * dV * dV;
  } // Close loop over cell m
  //  Energy = cahnHilliard + volumeConst + repulsion;
  cahnHilliard *= model->cahnHilliardCoeff;
  volumeConstraint *= model->volumePenaltyCoeff;
  repulsion = 0.5 * model->repulsionCoeff * (repulsion - phi4Sum);
  adhesion = 0.5 * model->adhesionCoeff * thickness2 * (adhesion - gphi2Sum);
  double energy = cahnHilliard + volumeConstraint + repulsion + adhesion;
  fprintf(f, "%d %g %g %g %g %g\n", step, cahnHilliard, volumeConstraint,
	  repulsion, adhesion, energy);
  fclose(f);
}

void deleteEnergyDump(EnergyDump* dump) {
  free(dump->phi2Field);
  free(dump->gphiField);
  free(dump);
}

DumpFuncs energyDumpFuncs =
  {
   .output = (void (*)(Dump*, Model*, int)) &energyOutput,
   .destroy = (void (*)(Dump*)) &deleteEnergyDump
  };

Dump* createEnergyDump(char* filename, int lx, int ly, int printInc) {
  EnergyDump* dump = malloc(sizeof *dump);
  setDump(&dump->super, dump, filename, printInc, &energyDumpFuncs);
  dump->lx = lx;
  dump->ly = ly;
  dump->phi2Field = create2DDoubleArray(lx, ly);
  dump->gphiField = create3DDoubleArray(lx, ly, 2);
  return (Dump*) dump;
}
