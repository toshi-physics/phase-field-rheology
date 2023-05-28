// model.h

#ifndef MODEL_H
#define MODEL_H

#include "cell.h"
#include "dump.h"

typedef struct Model {
  int lx;
  int ly;
  int cellLx;
  int cellLy;
  int numOfCells;
  int nshear;
  int doShear;
  double dt;
  double cellRadius;
  double thickness;
  double cahnHilliardCoeff;
  double volumePenaltyCoeff;
  double repulsionCoeff;
  double adhesionCoeff;
  double frictionCoeff;
  double diffusionCoeff;
  double activeShearCoeff;
  double mobility;
  double motility;
  double shearrate;
  double* cellXCM;
  double* cellYCM;
  int* cellXBoundCount;
  int* cellYBoundCount;
  Cell** cells;

  // Reduction fields
  double** totalField;
  double** totalField2;
  double** laplaceTotalField;
  double*** gradTotalField;
  double*** totalCellForceField;

  // Viscous terms
  double firstViscousCoeff;
  double secondViscousCoeff;
  
  // Dumps
  Dump** dumps;
  int ndumps;
  
  // Active and passive forces
  double* polarForce;
  double* totalForce; // Total force

  // Overlaps
  int doOverlap;
  double* viscousMat; // Overlap viscous matrix (in column major format)

  // Velocity from matrix inversion
  double* solvedVelocity; // Solved velocity vector for each cell
} Model;

Model* createModel(int lx, int ly,int ncells);
void deleteModel(Model* model);

void initCellsFromFile(Model* model, char* cellFile,
		       unsigned long seed);

void run(Model* model, int nsteps);
void output(Model* model, int step);

void updateTotalField(Model* model);
void updateTotalCellForceField(Model* model);
void updateViscous(Model* model);
void updateVelocity(Model* model);
void updateVelocityShear(Model* model);
void updateCellFieldGradient(Model* model, int m);
void updateCellChemPotAndDeform(Model* model, int m);
void updateCellPolarity(Model* model, int m);
void updateCellForces(Model* model, int m);
void updateCellField(Model* model, int m);
void updateCellCM(Model* model, int m);

#endif
