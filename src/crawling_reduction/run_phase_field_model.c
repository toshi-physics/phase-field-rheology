// run_phase_field.c

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "dump.h"
#include "phase_field_model.h"
#include "constant.h"

int main (int argc, char* argv[]) {
  
  if (argc != 2) {
    printf("Usage: RunPhaseFieldModel param_file [dump_file]\n");
    return 1;
  }

  int argi = 0;
  char* filename = argv[++argi];
  
  printf("Reading parameters ...\n");

  FILE* file = fopen(filename, "r");
  if (file == NULL) {
    printf("Error in opening parameter file!\n");
    return 1;
  }

  char line [PF_DIR_SIZE], dumpMode [PF_DIR_SIZE];
  char cmFile [PF_DIR_SIZE], shapeFile [PF_DIR_SIZE], dumpFile [PF_DIR_SIZE];
  double phi0 = -1.0;
  double mobility, cellRadius, dt, motility;
  double surfaceTensionCoeff, cahnHilliardCoeff, repulsionCoeff, adhesionCoeff,
    volumePenaltyCoeff, frictionCoeff, activeShearCoeff, diffusionCoeff;
  int cellLx = -1;
  int cellLy = -1;
  int lx = -1;
  int ly = -1;
  int nequil, nsteps, ncells;
  unsigned long seed;
  int nparams = 0;
  int ndumps = 0;
  int nedumps = 0;
  int maxDumps = 1000;
  int printInc, overwrite, cellIndex;
  int overlap;
  
#if PF_HAS_ARMA
  int fieldScale, kernelLength, sgolayDegree, sgolayLength;
  double kernelSigma;
#endif
  Dump** equilDumps = malloc(sizeof *equilDumps * maxDumps);
  Dump** dumps = malloc(sizeof *dumps * maxDumps);

  while (fgets(line, sizeof(line), file) != NULL) {
    nparams += sscanf(line, "phi0 = %lf", &phi0);
    nparams += sscanf(line, "cellRadius = %lf", &cellRadius);
    nparams += sscanf(line, "cahnHilliardCoeff = %lf", &cahnHilliardCoeff);
    nparams += sscanf(line, "surfaceTensionCoeff = %lf", &surfaceTensionCoeff);
    nparams += sscanf(line, "volumePenaltyCoeff = %lf", &volumePenaltyCoeff);
    nparams += sscanf(line, "repulsionCoeff = %lf", &repulsionCoeff);
    nparams += sscanf(line, "adhesionCoeff = %lf", &adhesionCoeff);
    nparams += sscanf(line, "frictionCoeff = %lf", &frictionCoeff);
    nparams += sscanf(line, "activeShearCoeff = %lf", &activeShearCoeff);
    nparams += sscanf(line, "diffusionCoeff = %lf", &diffusionCoeff);
    nparams += sscanf(line, "cellLx = %d", &cellLx);
    nparams += sscanf(line, "cellLy = %d", &cellLy);
    nparams += sscanf(line, "lx = %d", &lx);
    nparams += sscanf(line, "ly = %d", &ly);
    nparams += sscanf(line, "nsteps = %d", &nsteps);
    nparams += sscanf(line, "nequil = %d", &nequil);
    nparams += sscanf(line, "ncells = %d", &ncells);
    nparams += sscanf(line, "dt = %lf", &dt);
    nparams += sscanf(line, "mobility = %lf", &mobility);
    nparams += sscanf(line, "motility = %lf", &motility);
    nparams += sscanf(line, "cm_file = %s", cmFile);
    nparams += sscanf(line, "shape_file = %s", shapeFile);
    nparams += sscanf(line, "seed = %ld", &seed);
    nparams += sscanf(line, "overlap = %d", &overlap);
    
    // Read dumps
    // CM dump
    if (sscanf(line, "dump_cm %d %d %s %s", 
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	equilDumps[nedumps] = createCMDump(dumpFile, printInc, overwrite);
	nedumps++;
      } else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	dumps[ndumps] = createCMDump(dumpFile, printInc, overwrite);
	ndumps++; 
      }
    }
    // Bulk CM dump
    if (sscanf(line, "dump_bulk_cm %d %s %s",
	       &printInc, dumpMode, dumpFile) == 3) {
      if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	equilDumps[nedumps] = createBulkCMDump(dumpFile, printInc);
	nedumps++;
      } else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	dumps[ndumps] = createBulkCMDump(dumpFile, printInc);
	ndumps++; 
      }
    }
    // Gyration dump
    if (sscanf(line, "dump_gyr %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	equilDumps[nedumps] = 
	  createGyrationDump(dumpFile, printInc, overwrite);
	nedumps++;
      } else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	dumps[ndumps] = createGyrationDump(dumpFile, printInc, overwrite);
	ndumps++;
      }
    }
    // Gyration field dump
    if (sscanf(line, "dump_gyr_field %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	equilDumps[nedumps] = 
	  createGyrationFieldDump(dumpFile, printInc, overwrite);
	nedumps++;
      } else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	dumps[ndumps] = createGyrationFieldDump(dumpFile, printInc, overwrite);
	ndumps++;
      }
    }
    // Velocity dump
    if (sscanf(line, "dump_vel %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	equilDumps[nedumps] = 
	  createVelocityDump(dumpFile, printInc, overwrite);
	nedumps++;
      } else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	dumps[ndumps] = createVelocityDump(dumpFile, printInc, overwrite);
	ndumps++;
      }
    }
    // Velocity field dump
    if (sscanf(line, "dump_vel_field %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	equilDumps[nedumps] = 
	  createVelocityFieldDump(dumpFile, printInc, overwrite);
	nedumps++;
      } else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	dumps[ndumps] = createVelocityFieldDump(dumpFile, printInc, overwrite);
	ndumps++;
      }
    }
    // Deform dump
    if (sscanf(line, "dump_deform %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	equilDumps[nedumps] = 
	  createDeformDump(dumpFile, printInc, overwrite);
	nedumps++;
      } else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	dumps[ndumps] = createDeformDump(dumpFile, printInc, overwrite);
	ndumps++;
      }
    }
    // Deform field dump
    if (sscanf(line, "dump_deform_field %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	equilDumps[nedumps] = 
	  createDeformFieldDump(dumpFile, printInc, overwrite);
	nedumps++;
      } else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	dumps[ndumps] = createDeformFieldDump(dumpFile, printInc, overwrite);
	ndumps++;
      }
    }
    // Total field dump
    if (sscanf(line, "dump_field %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	equilDumps[nedumps] = createFieldDump(dumpFile, printInc, overwrite);
	nedumps++;
      } else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	dumps[ndumps] = createFieldDump(dumpFile, printInc, overwrite);
	ndumps++;
      }
    }
    // Total index field dump
    if (sscanf(line, "dump_index_field %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	equilDumps[nedumps] =
	  createIndexFieldDump(dumpFile, printInc, overwrite);
	nedumps++;
      } else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	dumps[ndumps] = createIndexFieldDump(dumpFile, printInc, overwrite);
	ndumps++;
      }
    }    
    // Individual cell field dump
    if (sscanf(line, "dump_cell_field %d %d %d %s %s",
	       &cellIndex, &printInc, &overwrite, dumpMode, dumpFile) == 5) {
      if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	equilDumps[nedumps] =
	  createCellFieldDump(dumpFile, cellIndex, printInc, overwrite);
	nedumps++;
      } else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	dumps[ndumps] =
	  createCellFieldDump(dumpFile, cellIndex, printInc, overwrite);
	ndumps++;
      }
    }
    // Neighbour dump
    if (sscanf(line, "dump_neighbour %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      // Only created the dump when the field size is known,
      // as it is needed for creating the neighbour analysers
      if (lx > 0 && ly > 0) {
	if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	  equilDumps[nedumps] =
	    createNeighbourDump(dumpFile, lx, ly, printInc, overwrite);
	  nedumps++;
	} else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	  dumps[ndumps] =
	    createNeighbourDump(dumpFile, lx, ly, printInc, overwrite);
	  ndumps++;
	}
      }
    }
    // Energy dump
    if (sscanf(line, "dump_energy %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      // Only created the dump when the field size is known,
      // as it is needed for creating the neighbour analysers
      if (lx > 0 && ly > 0) {
	if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	  equilDumps[nedumps] =
	    createEnergyDump(dumpFile, lx, ly, printInc, overwrite);
	  nedumps++;
	} else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	  dumps[ndumps] =
	    createEnergyDump(dumpFile, lx, ly, printInc, overwrite);
	  ndumps++;
	}
      }
    }
    // Overlap dump
    if (sscanf(line, "dump_overlap %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      // Only created the dump when the field size is known,
      // as it is needed for creating the overlap analysers
      if (cellLx > 0 && cellLy > 0) {
	if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	  equilDumps[nedumps] =
	    createOverlapDump(dumpFile, cellLx, cellLy, printInc, overwrite);
	  nedumps++;
	} else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	  dumps[ndumps] =
	    createOverlapDump(dumpFile, cellLx, cellLy, printInc, overwrite);
	  ndumps++;
	}
      }
    }
    // Overlap field dump
    if (sscanf(line, "dump_overlap_field %d %d %d %s %s",
	       &cellIndex, &printInc, &overwrite, dumpMode, dumpFile) == 5) {
      // Only created the dump when the field size is known,
      // as it is needed for creating the overlap analysers
      if (cellLx > 0 && cellLy > 0) {
	if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	  equilDumps[nedumps] = 
	    createOverlapFieldDump(dumpFile, cellLx, cellLy, cellIndex,
				   printInc, overwrite);
	  nedumps++;
	} else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	  dumps[ndumps] =
	    createOverlapFieldDump(dumpFile, cellLx, cellLy, cellIndex,
				   printInc, overwrite);
	  ndumps++;
	}
      }
    }
#if PF_HAS_ARMA
    // Shape dump
    if (sscanf(line, "dump_shape %d %d %lf %d %d %d %d %s %s",
	       &fieldScale, &kernelLength, &kernelSigma,
	       &sgolayDegree, &sgolayLength,
	       &printInc, &overwrite, dumpMode, dumpFile) == 9) {
      // Only created the dump when cell field size and phi0 are known,
      // as they are needed for creating the shape analysers
      if (cellLx > 0 && cellLy > 0 && phi0 >= 0.0) {
	if (strcmp(dumpMode, "equil") == 0 && nedumps+1 < maxDumps) {
	  equilDumps[nedumps] =
	    createShapeDump(dumpFile, fieldScale, cellLx, cellLy,
			    kernelLength, kernelSigma, sgolayDegree,
			    sgolayLength, phi0/2.0, printInc, overwrite);
	  nedumps++;
	} else if (strcmp(dumpMode, "main") == 0 && ndumps+1 < maxDumps) {
	  dumps[ndumps] =
	    createShapeDump(dumpFile, fieldScale, cellLx, cellLy,
			    kernelLength, kernelSigma, sgolayDegree,
			    sgolayLength, phi0/2.0, printInc, overwrite);
	  ndumps++;
	}
      }
    }
#endif
  }
  equilDumps = realloc(equilDumps, sizeof *equilDumps * nedumps);
  dumps = realloc(dumps, sizeof *dumps * ndumps);

  fclose(file);
  
  if (nparams != 24) {
    printf("Not enough parameters specified!\n");
    return 1;
  }  
  printf("Read parameters:\n");
  printf("lx = %d\n", lx);
  printf("ly = %d\n", ly);
  printf("cellLx = %d\n", cellLx);
  printf("cellLy = %d\n", cellLy);
  printf("ncells = %d\n", ncells);
  printf("phi0 = %.5f\n", phi0);
  printf("volumePenaltyCoeff = %.5f\n", volumePenaltyCoeff);
  printf("cahnHilliardCoeff = %.5f\n", cahnHilliardCoeff);
  printf("surfaceTensionCoeff = %.5f\n", surfaceTensionCoeff);
  printf("repulsionCoeff = %.5f\n", repulsionCoeff);
  printf("adhesionCoeff = %.5f\n", adhesionCoeff);
  printf("frictionCoeff = %.5f\n", frictionCoeff);
  printf("activeShearCoeff = %.5f\n", activeShearCoeff);
  printf("diffusionCoeff = %.5f\n", diffusionCoeff);
  printf("dt = %.5f\n", dt);
  printf("motility = %.5f\n", motility);
  printf("mobility = %.5f\n", mobility);
  printf("cellRadius = %.5f\n", cellRadius);
  printf("seed = %ld\n", seed);
  printf("cm_file = %s\n", cmFile);
  printf("shape_file = %s\n", shapeFile);  
  printf("nsteps = %d\n", nsteps);
  printf("nequil = %d\n", nequil);
  printf("overlap = %d\n", overlap);

  printf("Initialising model ...\n");
  
  PhaseFieldModel* model = createModel(lx, ly, ncells);
  model->phi0 = phi0;
  model->mobility = mobility;
  model->piR2phi02 = PF_PI * cellRadius * cellRadius * phi0 * phi0;
  model->surfaceTensionCoeff = surfaceTensionCoeff;
  model->cahnHilliardCoeff = cahnHilliardCoeff;
  model->volumePenaltyCoeff = volumePenaltyCoeff;
  model->diffusionCoeff = diffusionCoeff;
  model->repulsionCoeff = repulsionCoeff;
  model->adhesionCoeff = adhesionCoeff;
  model->frictionCoeff = frictionCoeff;
  model->activeShearCoeff = activeShearCoeff;
  model->motility = 0.0;
  model->cellLx = cellLx;
  model->cellLy = cellLy;
  model->ndumps = nedumps;
  model->dumps = equilDumps;
  model->doOverlap = overlap;
  
  initCellsFromFile(model, cmFile, shapeFile, seed);
  
  printf("Done initialisation.\n");
  
  model->dt = dt;

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

  model->motility = motility;
  model->ndumps = ndumps;
  model->dumps = dumps;

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
