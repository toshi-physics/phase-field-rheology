// phase_field_model.c

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "array.h"
#include "cell.h"
#include "phase_field_model.h"
#include "dump.h"
#include "random.h"
#include "constant.h"
#include "linear.h"
#include "util.h"

PhaseFieldModel* createModel(int lx, int ly, int ncells) {
  PhaseFieldModel* model =  malloc(sizeof *model);
  model->lx = lx;
  model->ly = ly;
  model->numOfCells = ncells;
  model->dt = 0.01;
  model->phi0 = 1.0;
  model->piR2phi02 = PF_PI;
  model->surfaceTensionCoeff = 0.0;
  model->cahnHilliardCoeff = 0.0;
  model->volumePenaltyCoeff = 0.0;
  model->repulsionCoeff = 0.0;
  model->adhesionCoeff = 0.0;
  model->frictionCoeff = 1.0;
  model->diffusionCoeff = 0.0;
  model->activeShearCoeff = 0.0;
  model->mobility = 1.0;
  model->motility = 0.0;
  model->cellLx = 1;
  model->cellLy = 1;
  model->cells = malloc(sizeof *model->cells * ncells);
  model->totalField = create2DDoubleArray(model->lx, model->ly);
  model->totalFieldSq = create2DDoubleArray(model->lx, model->ly);
  model->laplaceTotalField = create2DDoubleArray(model->lx, model->ly);
  model->totalCapillField = create3DDoubleArray(model->lx, model->ly, 2);
  model->totalDeformField = create3DDoubleArray(model->lx, model->ly, 2);
  model->cellXCM = create1DDoubleArray(ncells);
  model->cellYCM = create1DDoubleArray(ncells);
  model->cellXBoundCount = create1DIntArray(ncells);
  model->cellYBoundCount = create1DIntArray(ncells);
  model->dumps = NULL;
  model->ndumps = 0;
  model->overlapMat = create1DDoubleArray(ncells*ncells*4);
  model->overlapViscousMat = create1DDoubleArray(ncells*ncells*4);
  model->capillForce = create2DDoubleArray(ncells, 2);
  model->activeForce = create2DDoubleArray(ncells, 2);
  model->polarForce = create1DDoubleArray(ncells*2);
  model->totalForce = create1DDoubleArray(ncells*2);
  model->solvedVelocity = create1DDoubleArray(ncells*2);
  model->doOverlap = 1; // Enable overlap calculation by default
  model->firstViscousCoeff = 0.0;
  model->secondViscousCoeff = 0.0;
  model->viscousTensor = create3DDoubleArray(ncells, ncells, 4);
  return model; 
}

void deleteModel(PhaseFieldModel* model) {
  for (int i = 0; i < model->numOfCells; i++) {
    if (model->cells[i] != NULL) {
      deleteCell(model->cells[i]);
    }
  }
  free(model->cells);
  free(model->totalField);
  free(model->totalFieldSq);
  free(model->laplaceTotalField);
  free(model->totalCapillField);
  free(model->totalDeformField);
  free(model->cellXCM);
  free(model->cellYCM);
  free(model->cellXBoundCount);
  free(model->cellYBoundCount);
  free(model->overlapMat);
  free(model->overlapViscousMat);
  free(model->capillForce);
  free(model->activeForce);
  free(model->polarForce);
  free(model->totalForce);
  free(model->solvedVelocity);
  free(model->viscousTensor);
  free(model);
}

void initSquareCell(PhaseFieldModel* model, int index, int x, int y,
		    int dx, int dy) {
  int clx = model->cellLx;
  int cly = model->cellLy;
  Cell* cell = createCell(x, y, clx, cly, model->diffusionCoeff, 
			  model->phi0/2.0, index);
  model->cells[index] = cell;
  int x0 = (clx-dx)/2;
  int y0 = (cly-dy)/2;
  initFieldSquare(cell, x0, y0, dx, dy, model->phi0);
}

void initCellsFromFile(PhaseFieldModel* model, char* cmFile,
		       char* shapeFile, unsigned long seed) {
  FILE* fcm = fopen(cmFile, "r");
  if (fcm == NULL) {
    printf("Problem with opening the centre of mass file!\n");
    return;
  }

  FILE* fshape = fopen(shapeFile, "r");
  if (fshape == NULL) {
    printf("Problem with opening the shape file!\n");
    return;
  }

  char line [80];
  int nvar, x, y;
  double val;
  int clx = model->cellLx;
  int cly = model->cellLy;

  // Allocate memory for the template field
  double** field = create2DDoubleArray(clx, cly);
  // Read the template field from the shape file
  while (fgets(line, sizeof(line), fshape) != NULL) {
    nvar = sscanf(line, "%d %d %lf", &x, &y, &val);
    if (nvar == 3 && x >= 0 && x < clx && y >= 0 && y < cly) {
      field[x][y] = val;
    }
  }
  
  int index;
  double xcm, ycm, vx, vy;
  int count = 0;
  Cell* cell;
  while (fgets(line, sizeof(line), fcm) != NULL) {
    nvar = sscanf(line, "%d %lf %lf %lf %lf", &index, &xcm, &ycm, &vx, &vy);
    if (nvar == 3 || nvar == 5) {
      x = (int) round(xcm-clx/2.0);
      y = (int) round(ycm-cly/2.0);
      cell = createCell(x, y, clx, cly, model->diffusionCoeff, 
			model->phi0/2.0, index+seed);
      if (nvar == 5) {
	cell->vx = vx;
	cell->vy = vy;
	cell->v = sqrt(vx*vx+vy*vy);
	if (cell->v > 0.0) {
	  cell->theta = atan2(-vy,-vx)+PF_PI;
	}
      }
      model->cells[index] = cell;
      initField(cell, field);
      count++;
    }
  }

  if (count != model->numOfCells) {
    printf("Not all cells initialised!\n");
  }
  free(field);
  fclose(fcm);
  fclose(fshape);
}

void run(PhaseFieldModel* model, int nsteps) {
  output(model, 0); // Output initial state of the model
  for (int n = 1; n <= nsteps; n++) {
    updateTotalField(model);

    if (model->firstViscousCoeff > 0.0 || model->secondViscousCoeff > 0.0) {
#pragma omp parallel for default(none) shared(model) schedule(static)
      for (int m = 0; m < model->numOfCells; m++) {
	Cell* cell = model->cells[m];
	updateCellFieldGradient(cell);
      }
    }

    updateOverlap(model);

#pragma omp parallel for default(none) shared(model) schedule(static)
    for (int m = 0; m < model->numOfCells; m++) {
      Cell* cell = model->cells[m];
      updateCellPolarity(model, cell, m);
      updateCellChemPotAndDeform(model, cell, m);
    }
    
    updateTotalCapillDeformField(model);
    
#pragma omp parallel for default(none) shared(model) schedule(static)
    for (int m = 0; m < model->numOfCells; m++) {
      Cell* cell = model->cells[m];
      updateCellForces(model, cell, m);
    }
    
    updateVelocity(model);

#pragma omp parallel for default(none) shared(model) schedule(static)
    for (int m = 0; m < model->numOfCells; m++) {
      Cell* cell = model->cells[m];
      updateCellField(model, cell);
      updateCellCM(model, cell, m);
    }
    output(model, n);
  }
}

void output(PhaseFieldModel* model, int step) {
  if (step % 1000 == 0) {
    printf("Step %d\n", step);
  }
  for (int i = 0; i < model->ndumps; i++) {
    dumpOutput(model->dumps[i], model, step);
  }
}

void updateTotalField(PhaseFieldModel* model) {
  // Reset all global fields
  // No need to reset laplaceTotalField as values there are not summed
#pragma omp parallel for default(none) shared(model) schedule(static)
  for (int i = 0; i < model->lx; i++) {
    for (int j = 0; j < model->ly; j++) {
      model->totalField[i][j] = 0.0;
      model->totalFieldSq[i][j] = 0.0;
      model->totalCapillField[i][j][0] = 0.0;
      model->totalCapillField[i][j][1] = 0.0;
      model->totalDeformField[i][j][0] = 0.0;
      model->totalDeformField[i][j][1] = 0.0;
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
  shared(model, cell, clx, cly, cx, cy, cellField)	\
  private(x, y, phi) schedule(static)
    for (int i = 0; i < clx; i++) {
      x = iwrap(model->lx, cx+i);
      for (int j = 0; j < cly; j++) {
	y = iwrap(model->ly, cy+j);
	phi = cellField[i][j];
	model->totalField[x][y] += phi;
	model->totalFieldSq[x][y] += phi*phi;
      }
    }
  }

  int iu, id, ju, jd;
  
#pragma omp parallel for default(none)				\
  shared(model) private(iu, id, ju, jd) schedule(static)
  for (int i = 0; i < model->lx; i++) {
    iu = iup(model->lx, i);
    id = idown(model->lx, i);
    for (int j = 0; j < model->ly; j++) {
      ju = iup(model->ly, j);
      jd = idown(model->ly, j);
      model->laplaceTotalField[i][j] = 
	laplacian(i, j, iu, id, ju, jd, model->totalField);
    }
  }
}

void updateTotalCapillDeformField(PhaseFieldModel* model) {
  Cell* cell;
  double phi; 
  int clx, cly, x, y, cx, cy;
  double** cellField;

  for (int m = 0; m < model->numOfCells; m++) {
    cell = model->cells[m];
    clx = cell->lx;
    cly = cell->ly;
    cx = cell->x;
    cy = cell->y;
    cellField = cell->field[cell->getIndex];

#pragma omp parallel for default(none)			\
  shared(model, cell, clx, cly, cx, cy, cellField)	\
  private(x, y, phi) schedule(static)
    for (int i = 0; i < clx; i++) {
      x = iwrap(model->lx, cx+i);
      for (int j = 0; j < cly; j++) {
	y = iwrap(model->ly, cy+j);
	phi = cellField[i][j];
	model->totalCapillField[x][y][0] -= phi*cell->gradChemPot[i][j][0];
	model->totalCapillField[x][y][1] -= phi*cell->gradChemPot[i][j][1];
	model->totalDeformField[x][y][0] += cell->divDeform[i][j][0];
	model->totalDeformField[x][y][1] += cell->divDeform[i][j][1];
      }
    }
  }
}

void updateOverlap(PhaseFieldModel* model) {
  Cell* cellm;
  int ncells = model->numOfCells;
  int twoncells = ncells*2;
  double** phim;
  int clxm, clym;
  int doViscous = (model->firstViscousCoeff > 0.0 || 
		   model->secondViscousCoeff > 0.0);
  if (model->doOverlap || doViscous) {
    Cell* celln;
    double** phin;
    int x, y, cxm, cym, cxn, cyn, dxmn, dymn;
    double*** visTensor = model->viscousTensor;
    double* olapMat = model->overlapMat;
    double* olapVisMat = model->overlapViscousMat;
    // Reset the overlap matrix
    // This is LU factorised after each inversion, so must be reset
#pragma omp parallel for default(none)				\
  shared(olapMat, olapVisMat, ncells, doViscous) schedule(static)
    for (int m = 0; m < ncells*ncells*4; m++) {
      olapMat[m] = 0.0;
      if (doViscous) olapVisMat[m] = 0.0;
    }
    
    // Compute the actual amount of overlap \int dx \phi_i \phi_j
#pragma omp parallel for default(none) shared(model, ncells, twoncells) \
  shared(doViscous, olapMat, olapVisMat, visTensor)			\
  private(cellm, celln, phim, phin, x, y, clxm, clym, cxm, cym, cxn, cyn) \
  private(dxmn, dymn) schedule(dynamic, 4)
    for (int m = 0; m < ncells; m++) {
      cellm = model->cells[m];
      phim = cellm->field[cellm->getIndex];
      clxm = cellm->lx;
      clym = cellm->ly;
      cxm = iwrap(model->lx, cellm->x);
      cym = iwrap(model->ly, cellm->y);
      int twom = 2*m;
      int twomp1 = 2*m+1;
      
      // Compute self overlap
      double olap = 0.0;
      for (int i = 0; i < clxm; i++) {
	for (int j = 0; j < clym; j++) {
	  olap += phim[i][j]*phim[i][j];
	}
      }
      cellm->volume = olap;
      olapMat[twom*twoncells+twom] = olap;
      olapMat[twomp1*twoncells+twomp1] = olap;
      if (doViscous) {
	olapVisMat[twom*twoncells+twom] = olap;
	olapVisMat[twomp1*twoncells+twomp1] = olap;
      }
      
      for (int n = 0; n < m; n++) {
	celln = model->cells[n];
	cxn = iwrap(model->lx, celln->x);
	cyn = iwrap(model->ly, celln->y);
	dxmn = idiff(model->lx, cxm, cxn);
	dymn = idiff(model->ly, cym, cyn);
	int twon = 2*n;
	int twonp1 = 2*n+1;
	// Only do calculation if there is overlap 
	// between the two cells' domains
	if (abs(dxmn) < clxm && abs(dymn) < clym) {
	  phin = celln->field[celln->getIndex];
	  int xmstart, xnstart, ymstart, ynstart, lenx, leny, xm, ym, xn, yn;
	  // This assumes that clxm = clxn and clym = clyn
	  if (dxmn > 0) {
	    xmstart = 0;
	    xnstart = dxmn;
	    lenx = clxm-dxmn;
	  } else {
	    xmstart = -dxmn;
	    xnstart = 0;
	    lenx = clxm+dxmn;
	  }
	  if (dymn > 0) {
	    ymstart = 0;
	    ynstart = dymn;
	    leny = clym-dymn;
	  } else {
	    ymstart = -dymn;
	    ynstart = 0;
	    leny = clym+dymn;
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
	  double*** gphim = cellm->gradField;
	  double*** gphin = celln->gradField;
	  
	  for (int i = 0; i < lenx; i++) {
	    xm = xmstart+i;
	    xn = xnstart+i;
	    x = iwrap(model->lx, cxm+xm);
	    for (int j = 0; j < leny; j++) {
	      ym = ymstart+j;
	      yn = ynstart+j;
	      y = iwrap(model->ly, cym+ym);
	      
	      double pm = phim[xm][ym];
	      double pn = phin[xn][yn];
	      if (pm > 0.0 || pn > 0.0) {
		// Update overlap term
		olap += pm * pn;
	      
		// Update viscous term
		if (doViscous) {
		  double gpmx = gphim[xm][ym][0];
		  double gpmy = gphim[xm][ym][1];
		  double gpnx = gphin[xn][yn][0];
		  double gpny = gphin[xn][yn][1];
		  double tpsq = model->totalField[x][y] * 
		    model->totalField[x][y];
		  vismnxx += (pn * gpmx * gpmx - pm * gpmx * gpnx) / tpsq;
		  vismnxy += (pn * gpmx * gpmy - pm * gpmx * gpny) / tpsq;
		  vismnyx += (pn * gpmy * gpmx - pm * gpmy * gpnx) / tpsq;
		  vismnyy += (pn * gpmy * gpmy - pm * gpmy * gpny) / tpsq;
		  
		  visnmxx += (pm * gpnx * gpnx - pn * gpnx * gpmx) / tpsq;
		  visnmxy += (pm * gpnx * gpny - pn * gpnx * gpmy) / tpsq;
		  visnmyx += (pm * gpny * gpnx - pn * gpny * gpmx) / tpsq;
		  visnmyy += (pm * gpny * gpny - pn * gpny * gpmy) / tpsq;
		}
	      }
	    }
	  }
	  
	  olapMat[twom*twoncells+twon] = olap;
	  olapMat[twomp1*twoncells+twonp1] = olap;
	  olapMat[twon*twoncells+twom] = olap;
	  olapMat[twonp1*twoncells+twomp1] = olap;
	  
	  if (doViscous) {
	    double eta = model->firstViscousCoeff;
	    double zeta = model->secondViscousCoeff;
	    double vismnTrace = eta * (vismnxx + vismnyy);
	    double visnmTrace = eta * (visnmxx + visnmyy);
	    
	    vismnxx = vismnTrace + zeta * vismnxx;
	    vismnyy = vismnTrace + zeta * vismnyy;
	    vismnxy = zeta * vismnxy;
	    vismnyx = zeta * vismnyx;
	  
	    visnmxx = visnmTrace + zeta * visnmxx;
	    visnmyy = visnmTrace + zeta * visnmyy;
	    visnmxy = zeta * visnmxy;
	    visnmyx = zeta * visnmyx;
	    
	    olapVisMat[twon*twoncells+twom] = -vismnxx + olap;
	    olapVisMat[twon*twoncells+twomp1] = -vismnyx;
	    olapVisMat[twonp1*twoncells+twom] = -vismnxy;
	    olapVisMat[twonp1*twoncells+twomp1] = -vismnyy + olap;
	    
	    olapVisMat[twom*twoncells+twon] = -visnmxx + olap;
	    olapVisMat[twom*twoncells+twonp1] = -visnmyx;
	    olapVisMat[twomp1*twoncells+twon] = -visnmxy;
	    olapVisMat[twomp1*twoncells+twonp1] = -visnmyy + olap;

	    visTensor[m][n][0] = vismnxx;
	    visTensor[m][n][1] = vismnxy;
	    visTensor[m][n][2] = vismnyx;
	    visTensor[m][n][3] = vismnyy;
	    
	    visTensor[n][m][0] = visnmxx;
	    visTensor[n][m][1] = visnmxy;
	    visTensor[n][m][2] = visnmyx;
	    visTensor[n][m][3] = visnmyy;
	  }
	} // Close if block for overlap
      } // Close loop over cell n
    } // Close loop over cell m
    if (doViscous) {
#pragma omp parallel for default(none) \
  shared(ncells, twoncells, olapVisMat, visTensor) schedule(static)
      for (int m = 0; m < ncells; m++) {
	double vismxx = 0.0;
	double vismxy = 0.0;
	double vismyx = 0.0;
	double vismyy = 0.0;
	int twom = 2*m;
	int twomp1 = 2*m+1;
	for (int n = 0; n < ncells; n++) {
	  vismxx += visTensor[m][n][0];
	  vismxy += visTensor[m][n][1];
	  vismyx += visTensor[m][n][2];
	  vismyy += visTensor[m][n][3];
	}
	olapVisMat[twom*twoncells+twom] += vismxx;
	olapVisMat[twom*twoncells+twomp1] += vismyx;
	olapVisMat[twomp1*twoncells+twom] += vismxy;
	olapVisMat[twomp1*twoncells+twomp1] += vismyy;
      }
    }
  } else { // Don't do overlap and no viscous stress
#pragma omp parallel for default(none) shared(model, ncells) \
  private(cellm, phim, clxm, clym) 
    for (int m = 0; m < ncells; m++) {
      cellm = model->cells[m];
      clxm = cellm->lx;
      clym = cellm->ly;
      phim = cellm->field[cellm->getIndex];
      double olap = 0.0;
      for (int i = 0; i < clxm; i++) {
	for (int j = 0; j < clym; j++) {
	  olap += phim[i][j]*phim[i][j];
	}
      }
      cellm->volume = olap;
    }
  }
}

void updateCellFieldGradient(Cell* cell) {
  updateGradient(cell);
}

void updateCellChemPotAndDeform(PhaseFieldModel* model, Cell* cell, int m) {
  int clx = cell->lx;
  int cly = cell->ly;
  int cx = cell->x;
  int cy = cell->y;
  int x, y; // Lab frame coordinates of a lattice element
  int iu, iuu, id, idd, ju, juu, jd, jdd;  // Nearest neighbours
  double** cellField = cell->field[cell->getIndex];
  double*** gphi = cell->gradField;
  double phi, g2phi;
  double cahnHilliard, surfaceTension, volumeConst, repulsion, adhesion;
  // Need to update overlap first
  double vol = cell->volume / model->piR2phi02;
  
  // Apply fixed (Dirichlet) boundary conditions (phi = 0 at boundaries)
  // i and j are coordinates in the cell's own reference frame
  for (int i = 2; i < clx-2; i++) {
    iu = iup(clx, i);
    iuu = iup(clx, iu);
    id = idown(clx, i);
    idd = idown(clx, id);
    x = iwrap(model->lx, cx+i);
    for (int j = 2; j < cly-2; j++) {
      ju = iup(cly, j);
      juu = iup(cly, ju);
      jd = idown(cly, j);
      jdd = idown(cly, jd);
      y = iwrap(model->ly, cy+j);
      phi = cellField[i][j];
      g2phi = laplacian(i, j, iu, id, ju, jd, cellField);

      // Store the divergence of the deformation tensor
      // This is a shortcut that only works in 2D!
      cell->divDeform[i][j][0] = gphi[i][j][0]*g2phi;
      cell->divDeform[i][j][1] = gphi[i][j][1]*g2phi;

      // Cahn-Hilliard term
      cahnHilliard = model->cahnHilliardCoeff * phi * (phi - model->phi0) * 
	(phi - 0.5 * model->phi0);
      
      // Surface tension term
      surfaceTension = -model->surfaceTensionCoeff * g2phi;

      // Volume term
      volumeConst = 4.0 * model->volumePenaltyCoeff * phi * (vol - 1.0);
      
      // Repulsion term
      repulsion = 2.0 * model->repulsionCoeff * phi *
	(model->totalFieldSq[x][y] - phi*phi);
      
      // Adhesion term
      adhesion = model->adhesionCoeff * 
	(g2phi - model->laplaceTotalField[x][y]);
      
      // Store chemical potential
      cell->chemPot[i][j] = cahnHilliard + surfaceTension + volumeConst + 
	repulsion + adhesion;
    }
  }
  
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

void updateCellPolarity(PhaseFieldModel* model, Cell* cell, int m) {
  // Update polarity due to active rotational diffusion
  cell->theta += sqrt(6.0 * cell->diffusionCoeff * model->dt) * 
    (randDouble(cell->random)*2.0-1.0);
  int mx = 2*m;
  int my = mx+1;
  model->polarForce[mx] = model->motility * cos(cell->theta);
  model->polarForce[my] = model->motility * sin(cell->theta);
}

void updateCellForces(PhaseFieldModel* model, Cell* cell, int m) {
  // Compute the capillary force
  int x, y;
  double phi;
  int clx = cell->lx;
  int cly = cell->ly;
  int cx = cell->x;
  int cy = cell->y;
  double** cellField = cell->field[cell->getIndex];
  double* cap = model->capillForce[m];
  double* act = model->activeForce[m];
  double* tot = model->totalForce;
  cap[0] = 0.0;
  cap[1] = 0.0;
  act[0] = 0.0;
  act[1] = 0.0;
  for (int i = 0; i < clx; i++) {
    x = iwrap(model->lx, cx+i);
    for (int j = 0; j < cly; j++) {
      y = iwrap(model->ly, cy+j);
      phi = cellField[i][j];
      cap[0] += phi*model->totalCapillField[x][y][0];
      cap[1] += phi*model->totalCapillField[x][y][1];
      if (model->activeShearCoeff > 0.0) {
	act[0] += phi*model->totalDeformField[x][y][0];
	act[1] += phi*model->totalDeformField[x][y][1];
      }
    }
  }
  // Total force
  int mx = 2*m;
  int my = mx+1;
  tot[mx] = (cap[0] + model->activeShearCoeff*act[0]) / model->frictionCoeff;
  tot[my] = (cap[1] + model->activeShearCoeff*act[1]) / model->frictionCoeff;
}

void updateVelocity(PhaseFieldModel* model) {
  // This relies cell volume and polarity being up-to-date,
  // so update them first!
  int ncells = model->numOfCells;
  int twoncells = ncells*2;
  Cell* cell;
  int mx, my;
  double cvx, cvy;
  int error;
  int doViscous = (model->firstViscousCoeff > 0.0 || 
		   model->secondViscousCoeff > 0.0);
  if (model->doOverlap || doViscous) {
    ammpbm(model->overlapMat, model->polarForce, model->totalForce, 
	   1.0, 1.0, twoncells, twoncells, 1);
    if (doViscous) {
      error = solver(model->overlapViscousMat, model->totalForce, 
		     model->solvedVelocity, twoncells, 1, 0);
    } else {
      error = solver(model->overlapMat, model->totalForce, 
		     model->solvedVelocity, twoncells, 1, 0);
    }
    if (error != 0) {
      printf("Error: solution cannot be found for the matrix equation\n");
    }
    for (int m = 0; m < ncells; m++) {
      mx = 2*m;
      my = mx+1;
      cell = model->cells[m];
      cvx = model->solvedVelocity[mx];
      cvy = model->solvedVelocity[my];
      cell->vx = cvx;
      cell->vy = cvy;
      cell->v = sqrt(cvx*cvx+cvy*cvy);
    }
  } else { // Do simple division by cell volume
    for (int m = 0; m < ncells; m++) {
      mx = 2*m;
      my = mx+1;
      cell = model->cells[m];
      cvx = model->polarForce[mx] + model->totalForce[mx] / cell->volume;
      cvy = model->polarForce[my] + model->totalForce[my] / cell->volume;
      cell->vx = cvx;
      cell->vy = cvy;
      cell->v = sqrt(cvx*cvx+cvy*cvy);
    }
  }
}

void updateCellField(PhaseFieldModel* model, Cell* cell) {
  startUpdateCellField(cell);
  int clx = cell->lx;
  int cly = cell->ly;
  int set = cell->setIndex;
  int get = cell->getIndex;
  int iuu, iu, id, idd, juu, ju, jd, jdd;
  double vx = cell->vx;
  double vy = cell->vy;
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
      advection = upwind(i, j, iuu, iu, id, idd, 0, vx, cellField) +
	upwind(i, j, juu, ju, jd, jdd, 1, vy, cellField);

      cell->field[set][i][j] = cell->field[get][i][j] - model->dt *
	(model->mobility * cell->chemPot[i][j] + advection);
    }
  }
  endUpdateCellField(cell);
}

void updateCellCM(PhaseFieldModel* model, Cell* cell, int cellIndex) {
  updateCM(cell);
  int ix, iy;
  double x, y, cx, cy;
  cx = cell->x;
  cy = cell->y;
  x = cx+cell->xcm;
  y = cy+cell->ycm;
  ix = (int) floor(x / model->lx);
  iy = (int) floor(y / model->ly);
  model->cellXCM[cellIndex] = x - ix * model->lx;
  model->cellYCM[cellIndex] = y - iy * model->ly;
  model->cellXBoundCount[cellIndex] = ix;
  model->cellYBoundCount[cellIndex] = iy;
}
