// phase_field_model.h

#ifndef PHASE_FIELD_MODEL_H
#define PHASE_FIELD_MODEL_H

#include "cell.h"
#include "dump.h"

typedef struct PhaseFieldModel {
  int lx;
  int ly;
  int numOfCells;
  double dt;
  double piR2phi02;
  double phi0;
  double surfaceTensionCoeff;
  double cahnHilliardCoeff;
  double volumePenaltyCoeff;
  double diffusionCoeff;
  double mobility;
  double repulsionCoeff;
  double motility;
  double activeShearCoeff;
  double adhesionCoeff;
  int cellLx;
  int cellLy;
  double* cellXCM;
  double* cellYCM;
  int* cellXBoundCount;
  int* cellYBoundCount;
  Cell** cells;

  // Reduction fields
  double** totalField;
  double** totalFieldSq;
  double** laplaceTotalField;
  double** totalCapillField[2];
  double** totalDeformField[2]; // Traceless, active part

  // Dumps
  Dump** dumps;
  int ndumps;

  int* celldx; // Difference between top coordinates of cells in x
  int* celldy; // Difference between top coordinates of cells in y
  int* overlapInt; // Indicator if there is overlap or not
  double** overlap; // Actual overlap matrix (symmetric)
  double** capillForce; // Capillary force
  double** activeForce; // Active deviatoric force
  double* polarForce; // Polarity unity vector (in column major format)
  double* totalForce; // Total force (in column major format)
  double* overlapMat; // Overlap matrix (in column major format)
  double* solvedVelocity; // Solved velocity vector for each cell
} PhaseFieldModel;

PhaseFieldModel* createModel(int lx, int ly,int ncells);
void deleteModel(PhaseFieldModel* model);

void initSquareCell(PhaseFieldModel* model, int index,
		    int x, int y, int dx, int dy);

void initCellsFromFile(PhaseFieldModel* model, char* cmFile, char* shapeFile,
		       unsigned long seed);

void run(PhaseFieldModel* model, int nsteps);
void output(PhaseFieldModel* model, int step);

void updateTotalField(PhaseFieldModel* model);
void updateTotalCapillDeformField(PhaseFieldModel* model);
void updateOverlap(PhaseFieldModel* model);
void updateVelocity(PhaseFieldModel* model);
void updateCellDeform(PhaseFieldModel* model, Cell* cell);
void updateCellChemPotAndDeform(PhaseFieldModel* model, Cell* cell, int m);
//void updateCellVolume(PhaseFieldModel* model, Cell* cell);
void updateCellPolarity(PhaseFieldModel* model, Cell* cell, int m);
void updateCellForces(PhaseFieldModel* model, Cell* cell, int m);
void updateCellField(PhaseFieldModel* model, Cell* cell);
void updateCellCM(PhaseFieldModel* model, Cell* cell, int cellIndex);

int getPairIndex(int m, int n);

#endif
