// msd.cpp
// A code that computes the mean square displacement from the trajectory file

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <omp.h>
#include "position.hpp"
#include "array.hpp"

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::string;

int main (int argc, char* argv[]) {
  
  if (argc != 11) {
    cout << "Usage: msd npoints lx ly startTime endTime "
	 << "timeInc endShiftTime posFile posBulkFile outFile" << endl;
    return 1;
  }
  
  int argi = 0;
  int npoints = stoi(string(argv[++argi]), nullptr, 10);
  int lx = stoi(string(argv[++argi]), nullptr, 10);
  int ly = stoi(string(argv[++argi]), nullptr, 10);
  long startTime = stoi(string(argv[++argi]), nullptr, 10);
  long endTime = stoi(string(argv[++argi]), nullptr, 10);
  long timeInc = stoi(string(argv[++argi]), nullptr, 10);
  long endShiftTime = stoi(string(argv[++argi]), nullptr, 10);
  string posFile (argv[++argi]);
  string posBulkFile (argv[++argi]);
  string outFile (argv[++argi]);
  
  PositionReader reader;
  if (!reader.open(posFile, npoints, lx, ly, timeInc)) {
    cout << "Error: cannot open the file " << posFile << endl;
    return 1;
  }

  // Read the position data
  int nbins = static_cast<int>((endTime-startTime)/timeInc)+1;
  double*** pos = create3DArray<double>(nbins, npoints, 2);
  long time;
  int ibin;
  cout << "Reading data ..." << endl;
  while (reader.nextFrame()) {
    time = reader.getTime();
    if (time < startTime) {
      continue;
    } else if (time > endTime) {
      break;
    } else if ((time-startTime) % timeInc == 0) {
      ibin = static_cast<int>((time-startTime)/timeInc);
      for (int i = 0; i < npoints; i++) {
	pos[ibin][i][0] = reader.getUnwrappedPosition(i, 0);
	pos[ibin][i][1] = reader.getUnwrappedPosition(i, 1);
      }
    }
  }
  cout << "Done reading data" << endl;
  reader.close();

  // Read bulk cm
  double** totCM = create2DArray<double>(nbins, 2);
  ifstream cmReader;
  cmReader.open(posBulkFile);
  if (!cmReader) {
    cout << "Error: cannot open the file " << posBulkFile << endl;
    return 1;
  }
  long t;
  double xcm, ycm;
  stringstream iss;
  string line;
  while (getline(cmReader, line)) {
    iss.clear();
    iss.str(line);
    iss >> t >> xcm >> ycm;
    if (t < startTime) {
      continue;
    } else if (t > endTime) {
      break;
    } else if ((t-startTime) % timeInc == 0) {
      ibin = static_cast<int>((t-startTime)/timeInc);
      totCM[ibin][0] = xcm;
      totCM[ibin][1] = ycm;
    }
  }
  cmReader.close();

  int endShiftBin = static_cast<int>((endShiftTime-startTime)/timeInc)+1;
  if (endShiftBin >= nbins) endShiftBin = nbins-1;

#ifdef _OPENMP
  int nthreads = omp_get_max_threads();
#else
  int nthreads = 1;
#endif
  double** r2avg = create2DArray<double>(nthreads, nbins);
  double** r4avg = create2DArray<double>(nthreads, nbins);
  long** count = create2DArray<long>(nthreads, nbins);

#pragma omp parallel default(none) \
  shared(nbins, npoints, pos, r2avg, r4avg, count, totCM, endShiftBin)	\
  private(ibin)
  {
#ifdef _OPENMP
    int id = omp_get_thread_num();
#else
    int id = 0;
#endif
    double dx, dy, r2, r4, dxcm, dycm;
#pragma omp for schedule(dynamic,10)
    for (int i = 0; i < endShiftBin; i++) {
      for (int j = i+1; j < nbins; j++) {
	ibin = j-i;
	dxcm = totCM[j][0]-totCM[i][0];
	dycm = totCM[j][1]-totCM[i][1];
	for (int k = 0; k < npoints; k++) {
	  dx = (pos[j][k][0]-pos[i][k][0])-dxcm;
	  dy = (pos[j][k][1]-pos[i][k][1])-dycm;
	  r2 = dx*dx + dy*dy;
	  r4 = r2*r2;
	  r2avg[id][ibin] += r2;
	  r4avg[id][ibin] += r4;
	}
	count[id][ibin] += npoints;
      }     
    }
  }

  // Combine results
  for (int i = 1; i < nthreads; i++) {
    for (int j = 0; j < nbins; j++) {
      r2avg[0][j] += r2avg[i][j];
      r4avg[0][j] += r4avg[i][j];
      count[0][j] += count[i][j];
    }
  }
  // Normalise results
  for (int i = 1; i < nbins; i++) {
    r2avg[0][i] /= static_cast<double>(count[0][i]);
    r4avg[0][i] /= static_cast<double>(count[0][i]);
  }

  // Output results
  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Error: cannot open the file " << outFile << endl;
    return 1;
  }
  
  for (int i = 0; i < nbins; i++) {
    writer << (i*timeInc) << " " << r2avg[0][i] << " " << r4avg[0][i] << "\n";
  }
  writer.close();
  
  // Clean up
  deleteArray(pos);
  deleteArray(totCM);
  deleteArray(r2avg);
  deleteArray(r4avg);
  deleteArray(count);
}
