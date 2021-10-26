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
  double frictionCoeff;
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
  double*** totalCapillField;
  double*** totalDeformField; // Traceless, active part

  // Dumps
  Dump** dumps;
  int ndumps;
  
  double** capillForce; // Capillary force
  double** activeForce; // Active deviatoric force
  double* polarForce; // Polarity unity vector (in column major format)
  double* totalForce; // Total force (in column major format)
  int doOverlap; // A switch to turn on or off calculation on overlap
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
void updateCellChemPotAndDeform(PhaseFieldModel* model, Cell* cell, int m);
void updateCellPolarity(PhaseFieldModel* model, Cell* cell, int m);
void updateCellForces(PhaseFieldModel* model, Cell* cell, int m);
void updateCellField(PhaseFieldModel* model, Cell* cell);
void updateCellCM(PhaseFieldModel* model, Cell* cell, int cellIndex);

#endif
