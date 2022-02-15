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
  model->totalCapillField = create3DDoubleArray(model->lx, model->ly, 2);
  model->totalDeformField = create3DDoubleArray(model->lx, model->ly, 2);

  // Viscous terms
  model->firstViscousCoeff = 0.0;
  model->secondViscousCoeff = 0.0;
  model->viscousTensor = create3DDoubleArray(ncells, ncells, 4);

  // Dumps
  model->dumps = NULL;
  model->ndumps = 0;

  // Active and passive forces
  model->capillForce = create2DDoubleArray(ncells, 2);
  model->activeForce = create2DDoubleArray(ncells, 2);
  model->polarForce = create1DDoubleArray(ncells*2);
  model->totalForce = create1DDoubleArray(ncells*2);

  // Overlaps
  model->doOverlap = 1; // Enable overlap calculation by default
  model->overlapMat = create1DDoubleArray(ncells*ncells*4);
  model->overlapViscousMat = create1DDoubleArray(ncells*ncells*4);
  model->overlapTensor = create2DDoubleArray(ncells, ncells);

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
  free(model->overlapTensor);
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

    if (model->activeShearCoeff > 0.0 || model->firstViscousCoeff > 0.0 || 
	model->secondViscousCoeff > 0.0) {
#pragma omp parallel for default(none) shared(model) schedule(static)
      for (int m = 0; m < model->numOfCells; m++) {
	updateCellFieldGradient(model, m);
      }
    }

    updateOverlap(model);

#pragma omp parallel for default(none) shared(model) schedule(static)
    for (int m = 0; m < model->numOfCells; m++) {
      updateCellPolarity(model, m);
      updateCellChemPotAndDeform(model, m);
    }
    
    updateTotalCapillDeformField(model);
    
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
  if (model->activeShearCoeff > 0.0) {
#pragma omp parallel for default(none) shared(model) schedule(static)
    for (int i = 0; i < model->lx; i++) {
      for (int j = 0; j < model->ly; j++) {
	model->totalField[i][j] = 0.0;
	model->totalField2[i][j] = 0.0;
	model->totalCapillField[i][j][0] = 0.0;
	model->totalCapillField[i][j][1] = 0.0;
	model->totalDeformField[i][j][0] = 0.0;
	model->totalDeformField[i][j][1] = 0.0;
      }
    }
  } else {
#pragma omp parallel for default(none) shared(model) schedule(static)
    for (int i = 0; i < model->lx; i++) {
      for (int j = 0; j < model->ly; j++) {
	model->totalField[i][j] = 0.0;
	model->totalField2[i][j] = 0.0;
	model->totalCapillField[i][j][0] = 0.0;
	model->totalCapillField[i][j][1] = 0.0;
      }
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
	model->totalField2[x][y] += phi * phi;
      }
    }
  }

  // Only compute the Laplacian of the total field if there is adhesion
  if (model->adhesionCoeff > 0.0) {
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
}

void updateTotalCapillDeformField(Model* model) {
  Cell* cell;
  double phi; 
  int clx, cly, x, y, cx, cy;
  double** cellField;

  if (model->activeShearCoeff > 0.0) {
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
	  model->totalCapillField[x][y][0] -= phi * cell->gradChemPot[i][j][0];
	  model->totalCapillField[x][y][1] -= phi * cell->gradChemPot[i][j][1];
	  model->totalDeformField[x][y][0] += cell->divDeform[i][j][0];
	  model->totalDeformField[x][y][1] += cell->divDeform[i][j][1];
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
      cellField = cell->field[cell->getIndex];
      
#pragma omp parallel for default(none)			\
  shared(model, cell, clx, cly, cx, cy, cellField)	\
  private(x, y, phi) schedule(static)
      for (int i = 0; i < clx; i++) {
	x = iwrap(model->lx, cx+i);
	for (int j = 0; j < cly; j++) {
	  y = iwrap(model->ly, cy+j);
	  phi = cellField[i][j];
	  model->totalCapillField[x][y][0] -= phi * cell->gradChemPot[i][j][0];
	  model->totalCapillField[x][y][1] -= phi * cell->gradChemPot[i][j][1];
	}
      }
    }
  }
}

void updateOverlap(Model* model) {
  Cell* cellm;
  int ncells = model->numOfCells;
  int twoncells = ncells * 2;
  double** phim;
  int clxm, clym;
  int doViscous = (model->firstViscousCoeff > 0.0 || 
		   model->secondViscousCoeff > 0.0);
  if (model->doOverlap || doViscous) {
    Cell* celln;
    double** phin;
    int x, y, cxm, cym, cxn, cyn, dxmn, dymn;
    double*** visTensor = model->viscousTensor;
    double** olapTensor = model->overlapTensor;
    double* olapMat = model->overlapMat;
    double* olapVisMat = model->overlapViscousMat;
    // Reset the overlap matrix
    // This is LU factorised after each inversion, so must be reset
    /*#pragma omp parallel for default(none)				\
  shared(olapMat, olapVisMat, ncells, doViscous) schedule(static)
    for (int m = 0; m < ncells*ncells*4; m++) {
      olapMat[m] = 0.0;
      if (doViscous) olapVisMat[m] = 0.0;
      }*/
    
    // Compute the actual amount of overlap \int dx \phi_i \phi_j
#pragma omp parallel for default(none) shared(model, ncells, twoncells) \
  shared(doViscous, olapMat, olapVisMat, olapTensor, visTensor)		\
  private(cellm, celln, phim, phin, x, y, clxm, clym, cxm, cym, cxn, cyn) \
  private(dxmn, dymn) schedule(dynamic, 4)
    for (int m = 0; m < ncells; m++) {
      cellm = model->cells[m];
      phim = cellm->field[cellm->getIndex];
      clxm = cellm->lx;
      clym = cellm->ly;
      cxm = iwrap(model->lx, cellm->x);
      cym = iwrap(model->ly, cellm->y);
      
      // Compute self overlap
      double olap = 0.0;
      for (int i = 0; i < clxm; i++) {
	for (int j = 0; j < clym; j++) {
	  olap += phim[i][j] * phim[i][j];
	}
      }
      cellm->volume = olap;
      olapTensor[m][m] = olap;
      
      for (int n = 0; n < m; n++) {
	celln = model->cells[n];
	cxn = iwrap(model->lx, celln->x);
	cyn = iwrap(model->ly, celln->y);
	dxmn = idiff(model->lx, cxm, cxn);
	dymn = idiff(model->ly, cym, cyn);
	// Only do calculation if there is overlap 
	// between the two cells' domains
	if (abs(dxmn) < clxm && abs(dymn) < clym) {
	  phin = celln->field[celln->getIndex];
	  int xmstart, xnstart, ymstart, ynstart, lenx, leny, xm, ym, xn, yn;
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
	  double*** gphim = cellm->gradField;
	  double*** gphin = celln->gradField;
	  
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
	      if (pm > 0.0 || pn > 0.0) {
		// Update overlap term
		olap += pm * pn;
	      
		// Update viscous term
		if (doViscous) {
		  double gpmx = gphim[xm][ym][0];
		  double gpmy = gphim[xm][ym][1];
		  double gpnx = gphin[xn][yn][0];
		  double gpny = gphin[xn][yn][1];
		  double tphi2 = model->totalField[x][y] * 
		    model->totalField[x][y];
		  double pgmnx = (pn * gpmx - pm * gpnx) / tphi2;
		  double pgmny = (pn * gpmy - pm * gpny) / tphi2;

		  vismnxx += gpmx * pgmnx;
		  vismnxy += gpmx * pgmny;
		  vismnyx += gpmy * pgmnx;
		  vismnyy += gpmy * pgmny;
		  
		  visnmxx -= gpnx * pgmnx;
		  visnmxy -= gpnx * pgmny;
		  visnmyx -= gpny * pgmnx;
		  visnmyy -= gpny * pgmny;
		} // Close if statement for doViscous
	      } // Close if statement for non-zero pm * pn
	    } // Close loop over j
	  } // Close loop over i
	  
	  olapTensor[m][n] = olap;
	  
	  if (doViscous) {
	    double r2 = model->cellRadius * model->cellRadius;
	    double eta = model->firstViscousCoeff * r2;
	    double zeta = model->secondViscousCoeff * r2;
	    double vismnTrace = eta * (vismnxx + vismnyy);
	    double visnmTrace = eta * (visnmxx + visnmyy);
	    
	    visTensor[m][n][0] = -vismnTrace - zeta * vismnxx;
	    visTensor[m][n][1] = -zeta * vismnxy;
	    visTensor[m][n][2] = -zeta * vismnyx;
	    visTensor[m][n][3] = -vismnTrace - zeta * vismnyy;
	  
	    visTensor[n][m][0] = -visnmTrace - zeta * visnmxx;
	    visTensor[n][m][1] = -zeta * visnmxy;
	    visTensor[n][m][2] = -zeta * visnmyx;
	    visTensor[n][m][3] = -visnmTrace - zeta * visnmyy;
	  }
	} // Close if block for overlap
      } // Close loop over cell n
    } // Close loop over cell m (and close parallel region)
    
    //printf("Done computing viscous and overlap tensor\n");
    // Compute diagonal elements of the viscous tensor
    if (doViscous) {
#pragma omp parallel for default(none) \
  shared(ncells, visTensor) schedule(static)
      for (int m = 0; m < ncells; m++) {
	double vismxx = 0.0;
	double vismxy = 0.0;
	double vismyx = 0.0;
	double vismyy = 0.0;
	visTensor[m][m][0] = 0.0;
	visTensor[m][m][1] = 0.0;
	visTensor[m][m][2] = 0.0;
	visTensor[m][m][3] = 0.0;
	for (int n = 0; n < ncells; n++) {
	  vismxx -= visTensor[m][n][0];
	  vismxy -= visTensor[m][n][1];
	  vismyx -= visTensor[m][n][2];
	  vismyy -= visTensor[m][n][3];
	}
	visTensor[m][m][0] = vismxx;
	visTensor[m][m][1] = vismxy;
	visTensor[m][m][2] = vismyx;
	visTensor[m][m][3] = vismyy;
      }
    }
    
    // Combine overlap and viscous tensor (convert them into column major)
#pragma omp parallel for default(none) shared(ncells, twoncells, doViscous) \
  shared(visTensor, olapTensor, olapMat, olapVisMat) schedule(static)
    for (int m = 0; m < ncells; m++) {
      int twom = 2 * m;
      int twomp1 = 2 * m + 1;
      for (int n = 0; n < ncells; n++) {
	int twon = 2 * n;
	int twomn = twon * twoncells;
	int twomp1n = (twon + 1) * twoncells;
	int mnxx = twomn + twom;
	int mnyx = twomn + twomp1;
	int mnxy = twomp1n + twom;
	int mnyy = twomp1n + twomp1;
	double olap = (m > n ? olapTensor[m][n] : olapTensor[n][m]);
	olapMat[mnxx] = olap;
	olapMat[mnxy] = 0.0;
	olapMat[mnyx] = 0.0;
	olapMat[mnyy] = olap;
	if (doViscous) {
	  olapVisMat[mnxx] = (visTensor[m][n][0] + olap);
	  olapVisMat[mnxy] = visTensor[m][n][1];
	  olapVisMat[mnyx] = visTensor[m][n][2];
	  olapVisMat[mnyy] = (visTensor[m][n][3] + olap);
	}
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
	  olap += phim[i][j] * phim[i][j];
	}
      }
      cellm->volume = olap;
    }
  }
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
  int iu, iuu, id, idd, ju, juu, jd, jdd;  // Nearest neighbours
  double** cellField = cell->field[cell->getIndex];
  double*** gphi = cell->gradField;
  double phi, g2phi, cahnHilliard, volumeConstraint, repulsion;
  double adhesion = 0.0;
  double thickness2 = model->thickness * model->thickness;
  // Need to update overlap first
  double scaledVolume = cell->volume / (PF_PI * cell->radius * cell->radius);
  
  // Apply fixed (Dirichlet) boundary conditions (phi = 0 at boundaries)
  // i and j are coordinates in the cell's own reference frame
  for (int i = 2; i < clx-2; i++) {
    iu = iup(clx, i);
    id = idown(clx, i);
    x = iwrap(model->lx, cx+i);
    for (int j = 2; j < cly-2; j++) {
      ju = iup(cly, j);
      jd = idown(cly, j);
      y = iwrap(model->ly, cy+j);
      phi = cellField[i][j];
      g2phi = laplacian(i, j, iu, id, ju, jd, cellField);

      // Store the divergence of the deformation tensor
      // This is a shortcut that only works in 2D!
      if (model->activeShearCoeff > 0.0) {
	cell->divDeform[i][j][0] = gphi[i][j][0] * g2phi;
	cell->divDeform[i][j][1] = gphi[i][j][1] * g2phi;
      }

      // Cahn-Hilliard term
      cahnHilliard = 2.0 * model->cahnHilliardCoeff *
	(2.0 * phi * (phi - 1.0) * (phi - 0.5) - thickness2 * g2phi);
      
      // Volume term
      volumeConstraint = 4.0 * model->volumePenaltyCoeff * phi * 
	(scaledVolume - 1.0);
      
      // Repulsion term
      repulsion = 2.0 * model->repulsionCoeff * phi *
	(model->totalField2[x][y] - phi * phi);
      
      // Adhesion term
      if (model->adhesionCoeff > 0.0) {
	adhesion = model->adhesionCoeff * thickness2 *  
	  (g2phi - model->laplaceTotalField[x][y]);
      }
      
      // Store chemical potential
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
  int mx = 2 * m;
  int my = mx + 1;
  model->polarForce[mx] = model->motility * cos(cell->theta);
  model->polarForce[my] = model->motility * sin(cell->theta);
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
  double* cap = model->capillForce[m];
  double* tot = model->totalForce;
  cap[0] = 0.0;
  cap[1] = 0.0;
  if (model->activeShearCoeff > 0.0) {
    double* act = model->activeForce[m];
    act[0] = 0.0;
    act[1] = 0.0;
    for (int i = 0; i < clx; i++) {
      x = iwrap(model->lx, cx+i);
      for (int j = 0; j < cly; j++) {
	y = iwrap(model->ly, cy+j);
	phi = cellField[i][j];
	cap[0] += phi * model->totalCapillField[x][y][0];
	cap[1] += phi * model->totalCapillField[x][y][1];
	act[0] += phi * model->totalDeformField[x][y][0];
	act[1] += phi * model->totalDeformField[x][y][1];
      }
    }
    // Sum of capilliary and active deviatoric deformation force
    tot[mx] = (cap[0] / model->frictionCoeff + 
	       model->activeShearCoeff * act[0]);
    tot[my] = (cap[1] / model->frictionCoeff + 
	       model->activeShearCoeff * act[1]);
  } else {
    for (int i = 0; i < clx; i++) {
      x = iwrap(model->lx, cx+i);
      for (int j = 0; j < cly; j++) {
	y = iwrap(model->ly, cy+j);
	phi = cellField[i][j];
	cap[0] += phi * model->totalCapillField[x][y][0];
	cap[1] += phi * model->totalCapillField[x][y][1];
      }
    }
    tot[mx] = cap[0] / model->frictionCoeff;
    tot[my] = cap[1] / model->frictionCoeff;
  }
}

void updateVelocity(Model* model) {
  // This relies cell volume and polarity being up-to-date,
  // so update them first!
  int ncells = model->numOfCells;
  int twoncells = ncells * 2;
  Cell* cell;
  int mx, my;
  double cvx, cvy;
  int error;
  int doViscous = (model->firstViscousCoeff > 0.0 || 
		   model->secondViscousCoeff > 0.0);
  if (model->doOverlap || doViscous) {
    ammpbm(model->overlapMat, model->polarForce, model->totalForce, 
	   1.0, 1.0, twoncells, twoncells, 1);
    /*printMatrix(model->overlapMat, twoncells, twoncells);
    printf("\n");
    printMatrix(model->overlapViscousMat, twoncells, twoncells);
    printf("\n");*/
    if (doViscous) {
      error = solver(model->overlapViscousMat, model->totalForce, 
		     model->solvedVelocity, twoncells, 1, 0);
    } else {
      error = solver(model->overlapMat, model->totalForce, 
		     model->solvedVelocity, twoncells, 1, 0);
    }
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
  } else { // Do simple division by cell volume
    for (int m = 0; m < ncells; m++) {
      mx = 2 * m;
      my = mx + 1;
      cell = model->cells[m];
      cvx = model->polarForce[mx] + model->totalForce[mx] / cell->volume;
      cvy = model->polarForce[my] + model->totalForce[my] / cell->volume;
      cell->vx = cvx;
      cell->vy = cvy;
      cell->v = sqrt(cvx * cvx + cvy * cvy);
      //printf("%f ", cell->volume);
    }
    //printf("\n");
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
