// constant.h
// Stores the necessary constants in this program

#ifndef CONSTANT_H
#define CONSTANT_H

// Value for the mathemtical pi
#ifndef PF_PI
#define PF_PI 3.14159265358979323846264338327950288
#endif

// Width of halo on each side of cell subdomain
#ifndef PF_HALO
#define PF_HALO 2
#endif

// Threshold of change in cm before doing a frame shift
#ifndef PF_CMSHIFT
#define PF_CMSHIFT 2.0
#endif

#endif
