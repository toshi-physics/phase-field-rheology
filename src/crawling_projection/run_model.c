// run_model.c

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "dump.h"
#include "model.h"
#include "constant.h"

int main (int argc, char* argv[]) {
  
  if (argc != 2) {
    printf("Usage: run_model param_file\n");
    return 1;
  }

  int argi = 0;
  char* paramsFile = argv[++argi];
  
  printf("Reading parameters ...\n");

  FILE* fparams = fopen(paramsFile, "r");
  if (fparams == NULL) {
    printf("ERROR: cannot open the parameter file!\n");
    return 1;
  }

  char line [PF_DIR_SIZE], dumpMode [PF_DIR_SIZE];
  char cellFile [PF_DIR_SIZE], dumpFile [PF_DIR_SIZE];
  double cellRadius, thickness, cahnHilliardCoeff, volumePenaltyCoeff, 
    repulsionCoeff, adhesionCoeff, frictionCoeff, activeShearCoeff, 
    diffusionCoeff, firstViscousCoeff, secondViscousCoeff, dt, 
    mobility, motility;
  int nequil, nsteps, ncells, overlap;
  unsigned long seed;
  int lx = 0;
  int ly = 0;
  int cellLx = 0;
  int cellLy = 0;
  int nparams = 0;
  int ndumps = 0;
  int nedumps = 0;

  while (fgets(line, sizeof(line), fparams) != NULL) {
    nparams += sscanf(line, "cellRadius = %lf", &cellRadius);    
    nparams += sscanf(line, "thickness = %lf", &thickness);
    nparams += sscanf(line, "cahnHilliardCoeff = %lf", &cahnHilliardCoeff);
    nparams += sscanf(line, "volumePenaltyCoeff = %lf", &volumePenaltyCoeff);
    nparams += sscanf(line, "repulsionCoeff = %lf", &repulsionCoeff);
    nparams += sscanf(line, "adhesionCoeff = %lf", &adhesionCoeff);
    nparams += sscanf(line, "frictionCoeff = %lf", &frictionCoeff);
    nparams += sscanf(line, "activeShearCoeff = %lf", &activeShearCoeff);
    nparams += sscanf(line, "firstViscousCoeff = %lf", &firstViscousCoeff);
    nparams += sscanf(line, "secondViscousCoeff = %lf", &secondViscousCoeff);
    nparams += sscanf(line, "diffusionCoeff = %lf", &diffusionCoeff);
    nparams += sscanf(line, "lx = %d", &lx);
    nparams += sscanf(line, "ly = %d", &ly);
    nparams += sscanf(line, "cellLx = %d", &cellLx);
    nparams += sscanf(line, "cellLy = %d", &cellLy);
    nparams += sscanf(line, "nsteps = %d", &nsteps);
    nparams += sscanf(line, "nshear = %d", &nshear);
    nparams += sscanf(line, "nequil = %d", &nequil);
    nparams += sscanf(line, "ncells = %d", &ncells);
    nparams += sscanf(line, "dt = %lf", &dt);
    nparams += sscanf(line, "mobility = %lf", &mobility);
    nparams += sscanf(line, "motility = %lf", &motility);
    nparams += sscanf(line, "shearrate = %lf", &shearrate);
    nparams += sscanf(line, "cell_file = %s", cellFile);
    nparams += sscanf(line, "seed = %ld", &seed);
    nparams += sscanf(line, "overlap = %d", &overlap);
    
    // Count number of dumps
    if (line[0] != '#' && strstr(line, "dump_") != NULL) {
      if (strstr(line, " equil ") != NULL) {
	nedumps++;
      } else if (strstr(line, " main ") != NULL) {
	ndumps++;
      }
    }
  }

  fclose(fparams);

  if (nparams != 26) {
    printf("ERROR: Not enough parameters specified\n");
    return 1;
  }  
  printf("Read parameters:\n");
  printf("lx = %d\n", lx);
  printf("ly = %d\n", ly);
  printf("ncells = %d\n", ncells);
  printf("cellRadius = %g\n", cellRadius);
  printf("thickness = %g\n", thickness);
  printf("cahnHilliardCoeff = %g\n", cahnHilliardCoeff);
  printf("volumePenaltyCoeff = %g\n", volumePenaltyCoeff);
  printf("repulsionCoeff = %g\n", repulsionCoeff);
  printf("adhesionCoeff = %g\n", adhesionCoeff);
  printf("frictionCoeff = %g\n", frictionCoeff);
  printf("diffusionCoeff = %g\n", diffusionCoeff);
  printf("activeShearCoeff = %g\n", activeShearCoeff);
  printf("firstViscousCoeff = %g\n", firstViscousCoeff);
  printf("secondViscousCoeff = %g\n", secondViscousCoeff);
  printf("dt = %g\n", dt);
  printf("motility = %g\n", motility);
  printf("mobility = %g\n", mobility);
  printf("shearrate = %g\n", shearrate);
  printf("seed = %ld\n", seed);
  printf("overlap = %d\n", overlap);
  printf("cell_file = %s\n", cellFile);
  printf("nsteps = %d\n", nsteps);
  printf("nshear = %d\n", nshear);
  printf("nequil = %d\n", nequil);
  printf("nedumps = %d\n", nedumps);
  printf("ndumps = %d\n", ndumps);

  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  printf("nthreads = %d\n", nthreads);
  
  // Read dumps and fixes
  int idump = 0;
  int iedump = 0;
  int printInc, overwrite, cellIndex;
  int computeForceFreqEquil = -1;
  int computeForceFreq = -1;
  
#if PF_HAS_ARMA
  int fieldScale, kernelLength, sgolayDegree, sgolayLength, maxCellLx, 
    maxCellLy;
  double kernelSigma;
#endif
  Dump** equilDumps = malloc(sizeof *equilDumps * nedumps);
  Dump** dumps = malloc(sizeof *dumps * ndumps);

  fparams = fopen(paramsFile, "r");
  while (fgets(line, sizeof(line), fparams) != NULL) {
    // Read dumps
    // CM dump
    if (sscanf(line, "dump_cm %d %d %s %s", 
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0) {
	equilDumps[iedump] = createCMDump(dumpFile, printInc, overwrite);
	iedump++;
      } else if (strcmp(dumpMode, "main") == 0) {
	dumps[idump] = createCMDump(dumpFile, printInc, overwrite);
	idump++; 
      }
    }
    // Bulk CM dump
    if (sscanf(line, "dump_bulk_cm %d %s %s",
	       &printInc, dumpMode, dumpFile) == 3) {
      if (strcmp(dumpMode, "equil") == 0) {
	equilDumps[iedump] = createBulkCMDump(dumpFile, printInc);
	iedump++;
      } else if (strcmp(dumpMode, "main") == 0) {
	dumps[idump] = createBulkCMDump(dumpFile, printInc);
	idump++; 
      }
    }
    // Gyration dump
    if (sscanf(line, "dump_gyr %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0) {
	equilDumps[iedump] = 
	  createGyrationDump(dumpFile, printInc, overwrite);
	iedump++;
      } else if (strcmp(dumpMode, "main") == 0) {
	dumps[idump] = createGyrationDump(dumpFile, printInc, overwrite);
	idump++;
      }
    }
    // Gyration field dump
    if (sscanf(line, "dump_gyr_field %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0) {
	equilDumps[iedump] = 
	  createGyrationFieldDump(dumpFile, printInc, overwrite);
	iedump++;
      } else if (strcmp(dumpMode, "main") == 0) {
	dumps[idump] = createGyrationFieldDump(dumpFile, printInc, overwrite);
	idump++;
      }
    }
    // Velocity dump
    if (sscanf(line, "dump_vel %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0) {
	equilDumps[iedump] = 
	  createVelocityDump(dumpFile, printInc, overwrite);
	iedump++;
      } else if (strcmp(dumpMode, "main") == 0) {
	dumps[idump] = createVelocityDump(dumpFile, printInc, overwrite);
	idump++;
      }
    }
    // Velocity field dump
    if (sscanf(line, "dump_vel_field %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0) {
	equilDumps[iedump] = 
	  createVelocityFieldDump(dumpFile, printInc, overwrite);
	iedump++;
      } else if (strcmp(dumpMode, "main") == 0) {
	dumps[idump] = createVelocityFieldDump(dumpFile, printInc, overwrite);
	idump++;
      }
    }
    // Deform dump
    if (sscanf(line, "dump_deform %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0) {
	equilDumps[iedump] = 
	  createDeformDump(dumpFile, printInc, overwrite);
	iedump++;
      } else if (strcmp(dumpMode, "main") == 0) {
	dumps[idump] = createDeformDump(dumpFile, printInc, overwrite);
	idump++;
      }
    }
    // Deform field dump
    if (sscanf(line, "dump_deform_field %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0) {
	equilDumps[iedump] = 
	  createDeformFieldDump(dumpFile, printInc, overwrite);
	iedump++;
      } else if (strcmp(dumpMode, "main") == 0) {
	dumps[idump] = createDeformFieldDump(dumpFile, printInc, overwrite);
	idump++;
      }
    }
    // Total field dump
    if (sscanf(line, "dump_field %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0) {
	equilDumps[iedump] = createFieldDump(dumpFile, printInc, overwrite);
	iedump++;
      } else if (strcmp(dumpMode, "main") == 0) {
	dumps[idump] = createFieldDump(dumpFile, printInc, overwrite);
	idump++;
      }
    }
    // Total index field dump
    if (sscanf(line, "dump_index_field %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0) {
	equilDumps[iedump] =
	  createIndexFieldDump(dumpFile, printInc, overwrite);
	iedump++;
      } else if (strcmp(dumpMode, "main") == 0) {
	dumps[idump] = createIndexFieldDump(dumpFile, printInc, overwrite);
	idump++;
      }
    }    
    // Individual cell field dump
    if (sscanf(line, "dump_cell_field %d %d %d %s %s",
	       &cellIndex, &printInc, &overwrite, dumpMode, dumpFile) == 5) {
      if (strcmp(dumpMode, "equil") == 0) {
	equilDumps[iedump] =
	  createCellFieldDump(dumpFile, cellIndex, printInc, overwrite);
	iedump++;
      } else if (strcmp(dumpMode, "main") == 0) {
	dumps[idump] =
	  createCellFieldDump(dumpFile, cellIndex, printInc, overwrite);
	idump++;
      }
    }
    // Neighbour dump
    if (sscanf(line, "dump_neighbour %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      // Only created the dump when the field size is known,
      // as it is needed for creating the neighbour analysers
      if (lx > 0 && ly > 0) {
	if (strcmp(dumpMode, "equil") == 0) {
	  equilDumps[iedump] =
	    createNeighbourDump(dumpFile, lx, ly, printInc, overwrite);
	  iedump++;
	} else if (strcmp(dumpMode, "main") == 0) {
	  dumps[idump] =
	    createNeighbourDump(dumpFile, lx, ly, printInc, overwrite);
	  idump++;
	}
      }
    }
    // Force dump
    if (sscanf(line, "dump_force %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0) {
	equilDumps[iedump] = 
	  createForceDump(dumpFile, printInc, overwrite);
	computeForceFreqEquil = printInc;
	iedump++;
      } else if (strcmp(dumpMode, "main") == 0) {
	dumps[idump] = createForceDump(dumpFile, printInc, overwrite);
	computeForceFreq = printInc;
	idump++;
      }
    }
    // Energy dump
    if (sscanf(line, "dump_energy %d %s %s",
	       &printInc, dumpMode, dumpFile) == 3) {
      // Only created the dump when the field size is known,
      // as it is needed for creating the neighbour analysers
      if (lx > 0 && ly > 0) {
	if (strcmp(dumpMode, "equil") == 0) {
	  equilDumps[iedump] =
	    createEnergyDump(dumpFile, lx, ly, printInc);
	  iedump++;
	} else if (strcmp(dumpMode, "main") == 0) {
	  dumps[idump] =
	    createEnergyDump(dumpFile, lx, ly, printInc);
	  idump++;
	}
      }
    }
    // Overlap dump
    if (sscanf(line, "dump_overlap %d %d %d %d %s %s",
	       &printInc, &overwrite, &maxCellLx, &maxCellLy, 
	       dumpMode, dumpFile) == 6) {
      // Only created the dump when the field size is known,
      // as it is needed for creating the overlap analysers
      if (maxCellLx > 0 && maxCellLy > 0) {
	if (strcmp(dumpMode, "equil") == 0) {
	  equilDumps[iedump] = 
	    createOverlapDump(dumpFile, maxCellLx, maxCellLy, 
			      printInc, overwrite);
	  iedump++;
	} else if (strcmp(dumpMode, "main") == 0) {
	  dumps[idump] =
	    createOverlapDump(dumpFile, maxCellLx, maxCellLy, 
			      printInc, overwrite);
	  idump++;
	}
      }
    }
    // Overlap field dump
    if (sscanf(line, "dump_overlap_field %d %d %d %d %d %s %s",
	       &cellIndex, &printInc, &overwrite, &maxCellLx, &maxCellLy,
	       dumpMode, dumpFile) == 7) {
      // Only created the dump when the field size is known,
      // as it is needed for creating the overlap analysers
      if (maxCellLx > 0 && maxCellLy > 0) {
	if (strcmp(dumpMode, "equil") == 0) {
	  equilDumps[iedump] = 
	    createOverlapFieldDump(dumpFile, maxCellLx, maxCellLy, cellIndex,
				   printInc, overwrite);
	  iedump++;
	} else if (strcmp(dumpMode, "main") == 0) {
	  dumps[idump] =
	    createOverlapFieldDump(dumpFile, maxCellLx, maxCellLy, cellIndex,
				   printInc, overwrite);
	  idump++;
	}
      }
    }
#if PF_HAS_ARMA
    // Shape dump
    if (sscanf(line, "dump_shape %d %d %lf %d %d %d %d %d %d %s %s",
	       &fieldScale, &kernelLength, &kernelSigma,
	       &sgolayDegree, &sgolayLength, &maxCellLx, &maxCellLy, 
	       &printInc, &overwrite, dumpMode, dumpFile) == 11) {
      // Only created the dump when cell field size and phi0 are known,
      // as they are needed for creating the shape analysers
      if (maxCellLx > 0 && maxCellLy > 0) {
	if (strcmp(dumpMode, "equil") == 0) {
	  equilDumps[iedump] =
	    createShapeDump(dumpFile, fieldScale, maxCellLx, maxCellLy,
			    kernelLength, kernelSigma, sgolayDegree,
			    sgolayLength, 0.5, printInc, overwrite);
	  iedump++;
	} else if (strcmp(dumpMode, "main") == 0) {
	  dumps[idump] =
	    createShapeDump(dumpFile, fieldScale, maxCellLx, maxCellLy,
			    kernelLength, kernelSigma, sgolayDegree,
			    sgolayLength, 0.5, printInc, overwrite);
	  idump++;
	}
      }
    }
#endif
  } // Close loop over reading parameters
  
  fclose(fparams);
  
  printf("Initialising model ...\n");
  
  Model* model = createModel(lx, ly, ncells);
  model->cellLx = cellLx;
  model->cellLy = cellLy;
  model->cellRadius = cellRadius;
  model->thickness = thickness;
  model->cahnHilliardCoeff = cahnHilliardCoeff;
  model->volumePenaltyCoeff = volumePenaltyCoeff;
  model->repulsionCoeff = repulsionCoeff;
  model->adhesionCoeff = 0.0;
  model->frictionCoeff = 1.0;
  model->diffusionCoeff = diffusionCoeff;
  model->activeShearCoeff = 0.0;
  model->firstViscousCoeff = 0.0;
  model->secondViscousCoeff = 0.0;
  model->motility = 0.0;
  model->mobility = mobility;
  model->shearrate = shearrate;
  model->nshear = nshear;
  model->doShear = 0;
  model->ndumps = nedumps;
  model->dumps = equilDumps;
  model->dt = dt;
  model->doOverlap = overlap;
  model->computeForceFreq = computeForceFreqEquil;
  
  initCellsFromFile(model, cellFile, seed);
  
  printf("Done initialisation.\n");
  printf("Doing equilibration run ...\n");
  
#ifdef _OPENMP
  double start, end, duration;
  start = omp_get_wtime();
#endif
  run(model, nequil);

#ifdef _OPENMP
  end = omp_get_wtime();
  duration = end-start;
  printf("Time taken (sec): %.5f\n", duration);
  printf("\n");
#endif

  model->adhesionCoeff = adhesionCoeff;
  model->frictionCoeff = frictionCoeff;
  model->activeShearCoeff = activeShearCoeff;
  model->firstViscousCoeff = firstViscousCoeff;
  model->secondViscousCoeff = secondViscousCoeff;
  model->motility = motility;
  model->ndumps = ndumps;
  model->dumps = dumps;
  model->doShear = 1;
  model->computeForceFreq = computeForceFreq;
  
  printf("Doing main simulation run ...\n");

#ifdef _OPENMP
  start = omp_get_wtime();
#endif

  run(model, nsteps);

#ifdef _OPENMP
  end = omp_get_wtime();
  duration = end-start;
  printf("Time taken (sec): %.5f\n", duration);
#endif
  
  // Clean up
  deleteModel(model);

  for (int i = 0; i < nedumps; i++) {
    deleteDump(equilDumps[i]);
  }
  free(equilDumps);

  for (int i = 0; i < ndumps; i++) {
    deleteDump(dumps[i]);
  }
  free(dumps);
}
