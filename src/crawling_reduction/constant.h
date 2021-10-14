// constant.h
// Stores the necessary constants in this program

#ifndef CONSTANT_H
#define CONSTANT_H

// Value for the mathemtical pi
#ifndef PF_PI
#define PF_PI 3.14159265358979323846264338327950288
#endif

// Threshold of change in cm before doing a frame shift
#ifndef PF_CMSHIFT
#define PF_CMSHIFT 2.0
#endif

// Threhsold of overlap in phi_i^2phi_j^2 to be considered as neighbours
#ifndef PF_OLAPTHRES 
#define PF_OLAPTHRES 0.1
#endif

// Maximum number of neighbours per cell for neighbour list
#ifndef PF_MAXNEIGH
#define PF_MAXNEIGH 50
#endif

#endif
