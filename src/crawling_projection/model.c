// model.c

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "array.h"
#include "cell.h"
#include "model.h"
#include "dump.h"
#include "random.h"
#include "constant.h"
#include "linear.h"
#include "util.h"

Model* createModel(int lx, int ly, int ncells) {
  Model* model =  malloc(sizeof *model);
  model->lx = lx;
  model->ly = ly;
  model->cellLx = lx;
  model->cellLy = ly;
  model->numOfCells = ncells;
  model->dt = 0.01;
  model->cellRadius = 1.0;
  model->thickness = 1.0;
  model->cahnHilliardCoeff = 0.0;
  model->volumePenaltyCoeff = 0.0;
  model->repulsionCoeff = 0.0;
  model->adhesionCoeff = 0.0;
  model->frictionCoeff = 1.0;
  model->diffusionCoeff = 0.0;
  model->activeShearCoeff = 0.0;
  model->mobility = 1.0;
  model->motility = 0.0;
  model->cellXCM = create1DDoubleArray(ncells);
  model->cellYCM = create1DDoubleArray(ncells);
  model->cellXBoundCount = create1DIntArray(ncells);
  model->cellYBoundCount = create1DIntArray(ncells);
  model->cells = malloc(sizeof *model->cells * ncells);

  // Reduction fields
  model->totalField = create2DDoubleArray(model->lx, model->ly);
  model->totalField2 = create2DDoubleArray(model->lx, model->ly);
  model->laplaceTotalField = create2DDoubleArray(model->lx, model->ly);
  model->gradTotalField = create3DDoubleArray(model->lx, model->ly, 2);
  model->totalCellForceField = create3DDoubleArray(model->lx, model->ly, 2);

  // Viscous terms
  model->firstViscousCoeff = 0.0;
  model->secondViscousCoeff = 0.0;

  // Dumps
  model->dumps = NULL;
  model->ndumps = 0;

  // Force density
  model->totalForce = create1DDoubleArray(ncells*2);

  // Viscous matrix
  model->doOverlap  = 1;
  model->viscousMat = create1DDoubleArray(ncells*ncells*4);
  
  // Velocity from matrix inversion
  model->solvedVelocity = create1DDoubleArray(ncells*2);

  return model; 
}

void deleteModel(Model* model) {
  for (int i = 0; i < model->numOfCells; i++) {
    if (model->cells[i] != NULL) {
      deleteCell(model->cells[i]);
    }
  }
  free(model->cells);
  free(model->totalField);
  free(model->totalField2);
  free(model->laplaceTotalField);
  free(model->gradTotalField);
  free(model->totalCellForceField);
  free(model->cellXCM);
  free(model->cellYCM);
  free(model->cellXBoundCount);
  free(model->cellYBoundCount);
  free(model->viscousMat);
  free(model->totalForce);
  free(model->solvedVelocity);
  free(model);
}

void initCellsFromFile(Model* model, char* cellFile,
		       unsigned long seed) {
  FILE* fcell = fopen(cellFile, "r");
  if (fcell == NULL) {
    printf("Problem with opening the centre of mass file!\n");
    return;
  }
  
  char line [1000];
  int nvar, x, y, index;
  double xcm, ycm, vx, vy;
  int count = 0;
  int clx = model->cellLx;
  int cly = model->cellLy;
  Cell* cell;
  while (fgets(line, sizeof(line), fcell) != NULL) {
    nvar = sscanf(line, "%d %lf %lf %lf %lf", &index, &xcm, &ycm, &vx, &vy);
    if (nvar == 3 || nvar == 5) {
      x = (int) round(xcm-clx/2.0);
      y = (int) round(ycm-cly/2.0);
      cell = createCell(x, y, clx, cly, model->diffusionCoeff, 0.5, 
			index+seed);
      cell->radius = model->cellRadius;
      if (nvar == 5) {
	cell->vx = vx;
	cell->vy = vy;
	cell->v = sqrt(vx*vx+vy*vy);
	if (cell->v > 0.0) {
	  cell->theta = atan2(-vy,-vx)+PF_PI;
	  cell->px = cos(cell->theta);
	  cell->py = sin(cell->theta);
	}
      }
      model->cells[index] = cell;
      initCircleField(cell, model->cellRadius);
      count++;
    } else {
      printf("ERROR: not enough arguments supplied for a cell\n");
      exit(1);
    }
  }

  if (count != model->numOfCells) {
    printf("ERROR: Not all cells initialised\n");
    exit(1);
  }
  fclose(fcell);

  // Update all cells' CMs
  for (int m = 0; m < model->numOfCells; m++) {
    updateCellCM(model, m);
  }
}

void run(Model* model, int nsteps) {
  output(model, 0); // Output initial state of the model
  for (int n = 1; n <= nsteps; n++) {
    updateTotalField(model);

#pragma omp parallel for default(none) shared(model) schedule(static)
    for (int m = 0; m < model->numOfCells; m++) {
      updateCellFieldGradient(model, m);
    }
    
    updateViscous(model);

#pragma omp parallel for default(none) shared(model) schedule(static)
    for (int m = 0; m < model->numOfCells; m++) {
      updateCellPolarity(model, m);
      updateCellChemPotAndDeform(model, m);
    }
    
    updateTotalCellForceField(model);
    
#pragma omp parallel for default(none) shared(model) schedule(static)
    for (int m = 0; m < model->numOfCells; m++) {
      updateCellForces(model, m);
    }
    
    updateVelocity(model);

#pragma omp parallel for default(none) shared(model) schedule(static)
    for (int m = 0; m < model->numOfCells; m++) {
      updateCellField(model, m);
      updateCellCM(model, m);
    }
    output(model, n);
  }
}

void output(Model* model, int step) {
  if (step % 1000 == 0) {
    printf("Step %d\n", step);
  }
  for (int i = 0; i < model->ndumps; i++) {
    dumpOutput(model->dumps[i], model, step);
  }
}

void updateTotalField(Model* model) {
  // Reset all global fields
  // No need to reset laplaceTotalField as values there are not summed
  int lx = model->lx;
  int ly = model->ly;

#pragma omp parallel for default(none) shared(model, lx, ly) schedule(static)
  for (int i = 0; i < lx; i++) {
    for (int j = 0; j < ly; j++) {
      model->totalField[i][j] = 0.0;
      model->totalField2[i][j] = 0.0;
      model->totalCellForceField[i][j][0] = 0.0;
      model->totalCellForceField[i][j][1] = 0.0;
    }
  }

  Cell* cell;
  int clx, cly, x, y, cx, cy;
  double phi;
  double** cellField;
  
  for (int m = 0; m < model->numOfCells; m++) {
    cell = model->cells[m];
    clx = cell->lx;
    cly = cell->ly;
    cx = cell->x;
    cy = cell->y;
    cellField = cell->field[cell->getIndex];
#pragma omp parallel for default(none)			\
  shared(model, lx, ly, cell, clx, cly, cx, cy, cellField)	\
  private(x, y, phi) schedule(static)
    for (int i = 0; i < clx; i++) {
      x = iwrap(lx, cx+i);
      for (int j = 0; j < cly; j++) {
	y = iwrap(ly, cy+j);
	phi = cellField[i][j];
	model->totalField[x][y] += phi;
	model->totalField2[x][y] += phi * phi;
      }
    }
  }

  int iu, iuu, id, idd, ju, juu, jd, jdd;
  // Gradient of total field - needed for the viscous matrix
  if (model->firstViscousCoeff > 0.0 || model->secondViscousCoeff > 0.0) {
#pragma omp parallel for default(none) shared(model, lx, ly)	\
  private(iu, iuu, id, idd, ju, juu, jd, jdd) schedule(static)
    for (int i = 2; i < lx-2; i++) {
      iu = iup(lx, i);
      iuu = iup(lx, iu);
      id = idown(lx, i);
      idd = idown(lx, id);
      for (int j = 2; j < ly-2; j++) {
	ju = iup(ly, j);
	juu = iup(ly, ju);
	jd = idown(ly, j);
	jdd = idown(ly, jd);
	model->gradTotalField[i][j][0] = 
	  grad4(i, j, iuu, iu, id, idd, 0, model->totalField);
	model->gradTotalField[i][j][1] = 
	  grad4(i, j, juu, ju, jd, jdd, 1, model->totalField);
      }
    }
  }
  
  // Only compute the Laplacian of the total field if there is adhesion
  if (model->adhesionCoeff > 0.0) {
#pragma omp parallel for default(none) shared(model, lx, ly)	\
  private(iu, id, ju, jd) schedule(static)
    for (int i = 0; i < lx; i++) {
      iu = iup(lx, i);
      id = idown(lx, i);
      for (int j = 0; j < ly; j++) {
	ju = iup(ly, j);
	jd = idown(ly, j);
	model->laplaceTotalField[i][j] = 
	  laplacian(i, j, iu, id, ju, jd, model->totalField);
      }
    }
  }
}

void updateTotalCellForceField(Model* model) {
  Cell* cell;
  int clx, cly, x, y, cx, cy;
  double** cellField;
  double phi, tphi, mux, muy, px, py, cpx, cpy;
  double polarCoeff = model->frictionCoeff * model->motility;

  if (model->activeShearCoeff > 0.0) {
    double shearCoeff = model->activeShearCoeff * 
      model->thickness * model->thickness;
    double ax, ay;
    
    for (int m = 0; m < model->numOfCells; m++) {
      cell = model->cells[m];
      clx = cell->lx;
      cly = cell->ly;
      cx = cell->x;
      cy = cell->y;
      cpx = cell->px;
      cpy = cell->py;
      cellField = cell->field[cell->getIndex];
    
#pragma omp parallel for default(none) shared(model, cell, clx, cly, cx, cy) \
  shared(cpx, cpy, cellField, polarCoeff, shearCoeff)			\
  private(x, y, phi, tphi, mux, muy, px, py, ax, ay) schedule(static)
      for (int i = 0; i < clx; i++) {
	x = iwrap(model->lx, cx+i);
	for (int j = 0; j < cly; j++) {
	  y = iwrap(model->ly, cy+j);
	  tphi = model->totalField[x][y];
	  if (tphi > 0.0) {
	    phi = cellField[i][j];
	    mux = -phi * cell->gradChemPot[i][j][0];
	    muy = -phi * cell->gradChemPot[i][j][1];
	    px = phi * polarCoeff * cpx / tphi;
	    py = phi * polarCoeff * cpy / tphi;
	  }
	  ax = shearCoeff * cell->divDeform[i][j][0];
	  ay = shearCoeff * cell->divDeform[i][j][1];
	  model->totalCellForceField[x][y][0] += (mux + px + ax);
	  model->totalCellForceField[x][y][1] += (muy + py + ay);
	}
      }
    }
  } else {
    for (int m = 0; m < model->numOfCells; m++) {
      cell = model->cells[m];
      clx = cell->lx;
      cly = cell->ly;
      cx = cell->x;
      cy = cell->y;
      cpx = cell->px;
      cpy = cell->py;
      cellField = cell->field[cell->getIndex];

#pragma omp parallel for default(none)					\
  shared(model, cell, clx, cly, cx, cy, cpx, cpy, cellField, polarCoeff) \
  private(x, y, phi, tphi, mux, muy, px, py) schedule(static)
      for (int i = 0; i < clx; i++) {
	x = iwrap(model->lx, cx+i);
	for (int j = 0; j < cly; j++) {
	  y = iwrap(model->ly, cy+j);
	  tphi = model->totalField[x][y];
	  if (tphi > 0.0) {
	    phi = cellField[i][j];
	    mux = -phi * cell->gradChemPot[i][j][0];
	    muy = -phi * cell->gradChemPot[i][j][1];
	    px = phi * polarCoeff * cpx / tphi;
	    py = phi * polarCoeff * cpy / tphi;
	    model->totalCellForceField[x][y][0] += (mux + px);
	    model->totalCellForceField[x][y][1] += (muy + py);
	  }
	}
      }
    }
  }
}

void updateViscous(Model* model) {
  int ncells = model->numOfCells;
  int twoncells = ncells * 2;
  double eta = model->firstViscousCoeff;
  double zeta = model->secondViscousCoeff;
  double gamma = model->frictionCoeff;
  double* visMat = model->viscousMat;
  int doViscous = (eta > 0.0 || zeta > 0.0);

  // Compute the overlap 
#pragma omp parallel for default(none) shared(model, doViscous)		\
  shared(ncells, twoncells, visMat, eta, zeta, gamma) schedule(dynamic,1)
  for (int m = 0; m < ncells; m++) {
    Cell* cellm = model->cells[m];
    double** phim = cellm->field[cellm->getIndex];
    int clxm = cellm->lx;
    int clym = cellm->ly;
    int cxm = iwrap(model->lx, cellm->x);
    int cym = iwrap(model->ly, cellm->y);
    int twom = 2 * m;
    int twonm = twom * twoncells;

    // Compute self overlap
    double olap = 0.0;
    for (int i = 0; i < clxm; i++) {
      for (int j = 0; j < clym; j++) {
	olap += phim[i][j] * phim[i][j];
      }
    }
    cellm->volume = olap;
    
    if (model->doOverlap || doViscous) {
      for (int n = 0; n <= m; n++) {
	Cell* celln = model->cells[n];
	int cxn = iwrap(model->lx, celln->x);
	int cyn = iwrap(model->ly, celln->y);
	int dxmn = idiff(model->lx, cxm, cxn);
	int dymn = idiff(model->ly, cym, cyn);
	int twon = 2 * n;
	int twomn = twon * twoncells;
	int mnxx = twomn + twom;
	int mnyx = mnxx + 1;
	int mnxy = twomn + twoncells + twom;
	int mnyy = mnxy + 1;
	int nmxx = twonm + twon;
	int nmyx = nmxx + 1;
	int nmxy = twonm + twoncells + twon;
	int nmyy = nmxy + 1;
	  
	// Reset the matrix elements
	visMat[mnxx] = 0.0;
	visMat[mnyx] = 0.0;
	visMat[mnxy] = 0.0;
	visMat[mnyy] = 0.0;
	visMat[nmxx] = 0.0;
	visMat[nmyx] = 0.0;
	visMat[nmxy] = 0.0;
	visMat[nmyy] = 0.0;
	
	// Only do calculation if there is overlap 
	// between the two cells' domains
	if (abs(dxmn) < clxm && abs(dymn) < clym) {
	  int xmstart, xnstart, ymstart, ynstart, lenx, leny;
	  int xm, ym, xn, yn, x, y;	  
	  double** phin = celln->field[celln->getIndex];
	  // This assumes that clxm = clxn and clym = clyn
	  if (dxmn > 0) {
	    xmstart = 0;
	    xnstart = dxmn;
	    lenx = clxm - dxmn;
	  } else {
	    xmstart = -dxmn;
	    xnstart = 0;
	    lenx = clxm + dxmn;
	  }
	  if (dymn > 0) {
	    ymstart = 0;
	    ynstart = dymn;
	    leny = clym - dymn;
	  } else {
	    ymstart = -dymn;
	    ynstart = 0;
	    leny = clym + dymn;
	  }

	  olap = 0.0;

	  double vismnxx = 0.0;
	  double vismnyy = 0.0;
	  double vismnxy = 0.0;
	  double vismnyx = 0.0;
	  double visnmxx = 0.0;
	  double visnmyy = 0.0;
	  double visnmxy = 0.0;
	  double visnmyx = 0.0;
	  
	  if (doViscous) {
	    for (int i = 0; i < lenx; i++) {
	      xm = xmstart + i;
	      xn = xnstart + i;
	      x = iwrap(model->lx, cxm+xm);
	      for (int j = 0; j < leny; j++) {
		ym = ymstart + j;
		yn = ynstart + j;
		y = iwrap(model->ly, cym+ym);
		
		double pm = phim[xm][ym];
		double pn = phin[xn][yn];
		double tphi = model->totalField[x][y];
		
		if (tphi > 0.0) {
		  // Update overlap term
		  olap += pm * pn / tphi;
		  
		  // Update viscous term
		  double gpmx = cellm->gradField[xm][ym][0];
		  double gpmy = cellm->gradField[xm][ym][1];
		  double gpnx = celln->gradField[xn][yn][0];
		  double gpny = celln->gradField[xn][yn][1];
		  double gtpx = model->gradTotalField[x][y][0] / tphi;
		  double gtpy = model->gradTotalField[x][y][1] / tphi;
		  double vismx = (gpmx - pm * gtpx) / tphi;
		  double vismy = (gpmy - pm * gtpy) / tphi;
		  double visnx = (gpnx - pn * gtpx) / tphi;
		  double visny = (gpny - pn * gtpy) / tphi;
		  vismnxx += gpmx * visnx;
		  vismnxy += gpmx * visny;
		  vismnyx += gpmy * visnx;
		  vismnyy += gpmy * visny;	    
		  visnmxx += gpnx * vismx;
		  visnmxy += gpnx * vismy;
		  visnmyx += gpny * vismx;
		  visnmyy += gpny * vismy;
		}
	      } // Close loop over j
	    } // Close loop over i 
	    
	    olap *= gamma;
	    double vismndg = eta * (vismnxx + vismnyy) + olap;
	    double visnmdg = eta * (visnmxx + visnmyy) + olap;
	    
	    visMat[mnxx] = zeta * vismnxx + vismndg;
	    visMat[mnyx] = zeta * vismnyx;
	    visMat[mnxy] = zeta * vismnxy;
	    visMat[mnyy] = zeta * vismnyy + vismndg;
	    visMat[nmxx] = zeta * visnmxx + visnmdg;
	    visMat[nmyx] = zeta * visnmyx;
	    visMat[nmxy] = zeta * visnmxy;
	    visMat[nmyy] = zeta * visnmyy + visnmdg;
	  } else { // If eta and zeta are both zero
	    for (int i = 0; i < lenx; i++) {
	      xm = xmstart + i;
	      xn = xnstart + i;
	      x = iwrap(model->lx, cxm+xm);
	      for (int j = 0; j < leny; j++) {
		ym = ymstart + j;
		yn = ynstart + j;
		y = iwrap(model->ly, cym+ym);
		
		double pm = phim[xm][ym];
		double pn = phin[xn][yn];
		double tphi = model->totalField[x][y];
		// Update overlap term
		if (pm > 0.0 || pn > 0.0) {
		  olap += pm * pn / tphi;
		}
	      } // Close loop over j
	    } // Close loop over i
	    olap *= gamma;
	    visMat[mnxx] = olap;
	    visMat[mnyy] = olap;
	    visMat[nmxx] = olap;
	    visMat[nmyy] = olap;
	  } // Close if block for doViscous
	} // Close if block for overlap
      } // Close loop over cell n
    } // Close if block for doOverlap + doViscous
  } // Close loop over cell m (and close parallel region)
}

void updateCellFieldGradient(Model* model, int m) {
  updateGradient(model->cells[m]);
}

void updateCellChemPotAndDeform(Model* model, int m) {
  Cell* cell = model->cells[m];
  int clx = cell->lx;
  int cly = cell->ly;
  int cx = cell->x;
  int cy = cell->y;
  int x, y; // Lab frame coordinates of a lattice element
  int iu, iuu, id, idd, ju, juu, jd, jdd; // Nearest neighbours
  double** cellField = cell->field[cell->getIndex];
  double phi, lphi, cahnHilliard, volumeConstraint, repulsion;
  double adhesion = 0.0;
  double thickness2 = model->thickness * model->thickness;
  // Need to update overlap first
  double scaledVolume = cell->volume / (PF_PI * cell->radius * cell->radius);
  double athird = 1.0 / 3.0;
  double kappa = 2.0 * model->cahnHilliardCoeff;
  double epsilon = 2.0 * model->repulsionCoeff;
  double lambda = 4.0 * model->volumePenaltyCoeff;
  double omega = model->adhesionCoeff * thickness2;
  
  // Apply fixed (Dirichlet) boundary conditions (phi = 0 at boundaries)
  // i and j are coordinates in the cell's own reference frame
  for (int i = 2; i < clx-2; i++) {
    iu = iup(clx, i);
    id = idown(clx, i);
    if (model->activeShearCoeff > 0.0) {
      iuu = iup(clx, iu);
      idd = idown(clx, id);
    }
    x = iwrap(model->lx, cx+i);
    for (int j = 2; j < cly-2; j++) {
      ju = iup(cly, j);
      jd = idown(cly, j);
      y = iwrap(model->ly, cy+j);
      phi = cellField[i][j];
      lphi = laplacian(i, j, iu, id, ju, jd, cellField);
      
      // Compute divergence of deformation tensor
      if (model->activeShearCoeff > 0.0) {
	juu = iup(cly, ju);
	jdd = idown(cly, jd);
	double gpx = cell->gradField[i][j][0];
	double gpy = cell->gradField[i][j][1];
	double g2pxx = cgrad4(i, j, iuu, iu, id, idd, 0, 0, cell->gradField);
	double g2pxy = cgrad4(i, j, iuu, iu, id, idd, 1, 0, cell->gradField);
	double g2pyy = cgrad4(i, j, juu, ju, jd, jdd, 1, 1, cell->gradField);
	cell->divDeform[i][j][0] += 
	  athird * (g2pxx * gpx + g2pxy * gpy) + gpx * lphi;
	cell->divDeform[i][j][1] +=
	  athird * (g2pxy * gpx + g2pyy * gpy) + gpy * lphi; // g2pyx = g2pxy
      }

      // Compute chemical potential field
      // Cahn-Hilliard term
      cahnHilliard = kappa * 
	(2.0 * phi * (phi - 1.0) * (phi - 0.5) - thickness2 * lphi);
      
      // Volume term
      volumeConstraint = lambda * phi * (scaledVolume - 1.0);
      
      // Repulsion term
      repulsion = epsilon * phi * (model->totalField2[x][y] - phi * phi);
      
      // Adhesion term
      if (model->adhesionCoeff > 0.0) {
	adhesion = omega * (lphi - model->laplaceTotalField[x][y]);
      }
      
      cell->chemPot[i][j] = cahnHilliard + volumeConstraint + 
	repulsion + adhesion;
    } // Close loop over j
  } // Close loop over i
  
  // Compute gradient of the chemical potential
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
      cell->gradChemPot[i][j][0] = 
	grad4(i, j, iuu, iu, id, idd, 0, cell->chemPot);
      cell->gradChemPot[i][j][1] = 
	grad4(i, j, juu, ju, jd, jdd, 1, cell->chemPot);
    }
  }
}

void updateCellPolarity(Model* model, int m) {
  // Update polarity due to active rotational diffusion
  Cell* cell = model->cells[m];
  cell->theta += sqrt(6.0 * cell->diffusionCoeff * model->dt) * 
    (randDouble(cell->random) * 2.0 - 1.0);
  cell->px = cos(cell->theta);
  cell->py = sin(cell->theta);
}

void updateCellForces(Model* model, int m) {
  // Compute the capillary force and active deviatoric deformation force
  Cell* cell = model->cells[m];
  int x, y;
  double phi;
  int clx = cell->lx;
  int cly = cell->ly;
  int cx = cell->x;
  int cy = cell->y;
  int mx = 2 * m;
  int my = mx + 1;
  double** cellField = cell->field[cell->getIndex];
  double* tot = model->totalForce;
  double fx = 0.0;
  double fy = 0.0;
  for (int i = 0; i < clx; i++) {
    x = iwrap(model->lx, cx+i);
    for (int j = 0; j < cly; j++) {
      y = iwrap(model->ly, cy+j);
      phi = cellField[i][j];
      fx += phi * model->totalCellForceField[x][y][0];
      fy += phi * model->totalCellForceField[x][y][1];
    }
  }
  tot[mx] = fx;
  tot[my] = fy;
}

void updateVelocity(Model* model) {
  // This relies cell volume and polarity being up-to-date,
  // so update them first!
  int ncells = model->numOfCells;
  int twoncells = ncells * 2;
  Cell* cell;
  int mx, my;
  double cvx, cvy;
  int doViscous = (model->firstViscousCoeff > 0.0 || 
		   model->secondViscousCoeff > 0.0);
  if (model->doOverlap || doViscous) {
    int error = solver(model->viscousMat, model->totalForce, 
		       model->solvedVelocity, twoncells, 1, 0);
    if (error != 0) {
      printf("Error: solution cannot be found for the matrix equation\n");
      exit(1);
    }
    for (int m = 0; m < ncells; m++) {
      mx = 2 * m;
      my = mx + 1;
      cell = model->cells[m];
      cvx = model->solvedVelocity[mx];
      cvy = model->solvedVelocity[my];
      cell->vx = cvx;
      cell->vy = cvy;
      cell->v = sqrt(cvx * cvx + cvy * cvy);
    }
  } else { // Do simple divsion by cell volume
    for (int m = 0; m < ncells; m++) {
      mx = 2 * m;
      my = mx + 1;
      cell = model->cells[m];
      cvx = model->totalForce[mx] / cell->volume / model->frictionCoeff;
      cvy = model->totalForce[my] / cell->volume / model->frictionCoeff;
      cell->vx = cvx;
      cell->vy = cvy;
      cell->v = sqrt(cvx * cvx + cvy * cvy);
    }
  }
}

void updateCellField(Model* model, int m) {
  Cell* cell = model->cells[m];
  startUpdateCellField(cell);
  int clx = cell->lx;
  int cly = cell->ly;
  int set = cell->setIndex;
  int get = cell->getIndex;
  int iuu, iu, id, idd, juu, ju, jd, jdd;
  double advection;
  double** cellField = cell->field[get];
  
  // Apply fixed (Dirichlet) boundary conditions (phi = 0 at boundaries)
  // i and j are coordinates in the cell's own reference frame
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
      
      // Advection term (use the 3rd order upwind scheme)
      advection = upwind(i, j, iuu, iu, id, idd, 0, cell->vx, cellField) +
	upwind(i, j, juu, ju, jd, jdd, 1, cell->vy, cellField);

      cell->field[set][i][j] = cell->field[get][i][j] - model->dt *
	(model->mobility * cell->chemPot[i][j] + advection);
    }
  }
  endUpdateCellField(cell);
}

void updateCellCM(Model* model, int m) {
  Cell* cell = model->cells[m];
  updateCM(cell);
  int ix, iy;
  double x, y, cx, cy;
  cx = cell->x;
  cy = cell->y;
  x = cx + cell->xcm;
  y = cy + cell->ycm;
  ix = (int) floor(x / model->lx);
  iy = (int) floor(y / model->ly);
  model->cellXCM[m] = x - ix * model->lx;
  model->cellYCM[m] = y - iy * model->ly;
  model->cellXBoundCount[m] = ix;
  model->cellYBoundCount[m] = iy;
}
