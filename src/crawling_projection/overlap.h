// overlap.h
// Code to analyse the degree of overlap between phase fields

#ifndef OVERLAP_H
#define OVERLAP_H

#include "model.h"

typedef struct OverlapAnalyser {
  int numOfCells;
  int cellLx;
  int cellLy;
  double** overlapField;
} OverlapAnalyser; 

OverlapAnalyser* createOverlapAnalyser(int clx, int cly);
void deleteOverlapAnalyser(OverlapAnalyser* ana);
double computeOverlap(OverlapAnalyser* ana, Model* model, int cellIndex);

#endif
