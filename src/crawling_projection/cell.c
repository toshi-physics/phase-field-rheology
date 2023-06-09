// cell.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "cell.h"
#include "array.h"
#include "random.h"
#include "constant.h"
#include "util.h"

Cell* createCell(int x, int y, int lx, int ly,
		 double dr, double incell, unsigned long seed) {
  // Allocate memory for a cell
  Cell* cell = malloc(sizeof *cell);
  
  cell->x = x;
  cell->y = y;
  cell->lx = lx;
  cell->ly = ly;
  cell->incell = incell;

  // Allocate memory and initialise the field to zero
  for (int i = 0; i < 2; i++) {
    cell->field[i] = create2DDoubleArray(cell->lx, cell->ly);
  }
  cell->setIndex = 1;
  cell->getIndex = 1;
  cell->diffusionCoeff = dr;
  cell->xcm = cell->lx/2.0;
  cell->ycm = cell->ly/2.0;
  cell->deltaXCM = 0.0;
  cell->deltaYCM = 0.0;
  cell->volume = 0.0;
  cell->radius = 1.0;
  cell->random = createRandom(seed);
  // pick a random 2D direction
  cell->theta = randDouble(cell->random)*2.0*PF_PI;
  double costheta = cos(cell->theta);
  double sintheta = sin(cell->theta);
  cell->px = costheta;
  cell->py = sintheta;
  cell->vx = costheta;
  cell->vy = sintheta;
  cell->v = 1.0;
  cell->chemPot = create2DDoubleArray(cell->lx, cell->ly);
  cell->gradChemPot = create3DDoubleArray(cell->lx, cell->ly, 2);
  cell->divDeform = create3DDoubleArray(cell->lx, cell->ly, 2);
  cell->gradField = create3DDoubleArray(cell->lx, cell->ly, 2);
  return cell;
}

void deleteCell(Cell* cell) {
  for (int i = 0; i < 2; i++) {
    free(cell->field[i]);
  }
  deleteRandom(cell->random);
  free(cell->chemPot);
  free(cell->gradChemPot);
  free(cell->divDeform);
  free(cell->gradField);
  free(cell);
}

void initCircleField(Cell* cell, double r) {
  double x0 = cell->lx/2.0;
  double y0 = cell->ly/2.0;
  double r2 = r*r;
  double dx, dy;
  for (int i = 0; i < cell->lx; i++) {
    dx = i+0.5-x0;
    for (int j = 0; j < cell->ly; j++) {
      dy = j+0.5-y0;
      if (dx*dx+dy*dy <= r2) {
	cell->field[0][i][j] = 1.0;
	cell->field[1][i][j] = 1.0;
      }
    }
  }
  calculateCM(cell, &cell->xcm, &cell->ycm);
  updateVolume(cell);
}

void calculateCM(Cell* cell, double* xcm, double* ycm) {
  double xavg = 0.0;
  double yavg = 0.0;
  double mass = 0.0;
  double phi;
  int get = cell->getIndex;
  for (int i = 0; i < cell->lx; i++) {
    for (int j = 0; j < cell->ly; j++) {
      phi = cell->field[get][i][j];
      xavg += (phi * (i+0.5)); // Use the centre of a lattice element
      yavg += (phi * (j+0.5));
      mass += phi;
    }
  }
  if (mass > 0.0) {
    xavg /= mass;
    yavg /= mass;
  } else {
    xavg = 0.0;
    yavg = 0.0;
  }
  *xcm = xavg;
  *ycm = yavg;
}

void updateCM(Cell* cell) {
  double oldXCM = cell->xcm;
  double oldYCM = cell->ycm;
  calculateCM(cell, &cell->xcm, &cell->ycm);
  cell->drx = cell->xcm - oldXCM;
  cell->dry = cell->ycm - oldYCM;
  cell->deltaXCM += (cell->xcm - oldXCM);
  cell->deltaYCM += (cell->ycm - oldYCM);

  if (fabs(cell->deltaXCM) > PF_CMSHIFT || fabs(cell->deltaYCM) > PF_CMSHIFT) {
    int xShift = (int)(round(cell->deltaXCM));
    int yShift = (int)(round(cell->deltaYCM));
    shiftCoordinates(cell, xShift, yShift);
    cell->deltaXCM -= xShift;
    cell->deltaYCM -= yShift;
    cell->xcm -= xShift;
    cell->ycm -= yShift;
    cell->x += xShift;
    cell->y += yShift;    
  }
}

void updateVolume(Cell* cell) {
  double totalVolume = 0.0;
  int get = cell->getIndex;
  double phi;
  for (int i = 2; i < cell->lx-2; i++) {
    for (int j = 2; j < cell->ly-2; j++) {
      phi = cell->field[get][i][j];
      totalVolume += phi*phi;
    }
  }
  cell->volume = totalVolume;
}

void updateGradient(Cell* cell) {
  int iuu, iu, id, idd, juu, ju, jd, jdd;
  int clx = cell->lx;
  int cly = cell->ly;
  double** field = cell->field[cell->getIndex];
  for (int i = 0; i < cell->lx; i++) {
    iu = iup(clx, i);
    iuu = iup(clx, iu);
    id = idown(clx, i);
    idd = idown(clx, id);
    for (int j = 0; j < cell->ly; j++) {
      ju = iup(cly, j);
      juu = iup(cly, ju);
      jd = idown(cly, j);
      jdd = idown(cly, jd);
      cell->gradField[i][j][0] = grad4(i, j, iuu, iu, id, idd, 0, field);
      cell->gradField[i][j][1] = grad4(i, j, juu, ju, jd, jdd, 1, field);
    }
  }
}

void shiftCoordinates(Cell* cell, int xShift, int yShift) {
  startUpdateCellField(cell);
  int set = cell->setIndex;
  int get = cell->getIndex;
  // These are all in the new frame
  int clx = cell->lx;
  int cly = cell->ly;
  int xStart, xEnd, yStart, yEnd;
  int zeroXStart, zeroXEnd, zeroYStart, zeroYEnd;
  if (xShift >= 0) {
    xStart = 0;
    xEnd = clx - xShift;
    zeroXStart = xEnd;
    zeroXEnd = clx;
  } else {
    xStart = -xShift;
    xEnd = clx;
    zeroXStart = 0;
    zeroXEnd = xStart;
  }
  if (yShift >= 0) {
    yStart = 0;
    yEnd = cly - yShift;
    zeroYStart = yEnd;
    zeroYEnd = cly;
  } else {
    yStart = -yShift;
    yEnd = cly;
    zeroYStart = 0;
    zeroYEnd = yStart;
  }
  // Copy the data
  for (int i = xStart; i < xEnd; i++) {
    for (int j = yStart; j < yEnd; j++) {
      cell->field[set][i][j] = cell->field[get][i+xShift][j+yShift];
    }
  }
  // Set empty cells to zero
  for (int i = zeroXStart; i < zeroXEnd; i++) {
    for (int j = 0; j < cly; j++) {
      cell->field[set][i][j] = 0.0;
    }
  }
  if (zeroYStart < zeroYEnd) {
    for (int i = 0; i < clx; i++) {
      for (int j = zeroYStart; j < zeroYEnd; j++) {
	cell->field[set][i][j] = 0.0;
      }
    }
  }
  endUpdateCellField(cell);
}

inline void startUpdateCellField(Cell* cell) {
  cell->setIndex = (cell->getIndex == 1 ? 0 : 1);
}

inline void endUpdateCellField(Cell* cell) {
  cell->getIndex = cell->setIndex;
}
