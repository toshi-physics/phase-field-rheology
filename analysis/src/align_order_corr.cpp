// align_order_corr.cpp
// Output the pairwise alignment order and the radial distance

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include "position.hpp"
#include "array.hpp"

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::string;

// Helper function
double sgnd(double val);
double ddiff(double len, double d1, double d2);

int main(int argc, char* argv[]) {
  if (argc != 13) {
    cout << "Usage: align_order_corr npoints lx ly startTime endTime timeInc "
	 << "rmin rmax rinc posFile vecFile outFile" << endl;
    return 1;
  }

  int argi = 0;
  int npoints = stoi(string(argv[++argi]), nullptr, 10);
  int lx = stoi(string(argv[++argi]), nullptr, 10);
  int ly = stoi(string(argv[++argi]), nullptr, 10);
  long startTime = stol(string(argv[++argi]), nullptr, 10);
  long endTime = stol(string(argv[++argi]), nullptr, 10);
  long timeInc = stol(string(argv[++argi]), nullptr, 10);
  double rmin = stod(string(argv[++argi]), nullptr);
  double rmax = stod(string(argv[++argi]), nullptr);
  double rinc = stod(string(argv[++argi]), nullptr);
  string posFile (argv[++argi]);
  string vecFile (argv[++argi]);
  string outFile (argv[++argi]);

  PositionReader posReader;
  if (!posReader.open(posFile, npoints, lx, ly, timeInc)) {
    cout << "Error: cannot open the file " << posFile << endl;
    return 1;
  }
  
  // Read the position data
  int nbins = static_cast<int>((endTime-startTime)/timeInc)+1;
  double*** pos = create3DArray<double>(nbins, npoints, 2);
  long time;
  int ibin;
  cout << "Reading position data ..." << endl;
  while (posReader.nextFrame()) {
    time = posReader.getTime();
    if (time < startTime) {
      continue;
    } else if (time > endTime) {
      break;
    } else if ((time-startTime) % timeInc == 0) {
      ibin = static_cast<int>((time-startTime)/timeInc);
      for (int i = 0; i < npoints; i++) {
	pos[ibin][i][0] = posReader.getPosition(i, 0);
	pos[ibin][i][1] = posReader.getPosition(i, 1);
      }
    }
  }
  posReader.close();
  
  ifstream vecReader;
  vecReader.open(vecFile);
  if (!vecReader) {
    cout << "Error: cannot open the file " << vecFile << endl;
    return 1;
  }
  
  string line, str;
  stringstream ss;
  double** align = create2DArray<double>(npoints, 2);
  int rnbins = static_cast<int>(ceil((rmax-rmin)/rinc));
  double** alignCorr = create2DArray<double>(rnbins, 2);
  int* alignCount = create1DArray<int>(rnbins);
  int rbin;
  double w, vx, vy, dx, dy, dr, w2avg, cos2dt, wcos2dt;
  while (getline(vecReader, line)) {
    // Read the two header lines and get time
    getline(vecReader, line);
    ss.clear();
    ss.str(line);
    ss >> str >> time;
    // Only use the data from the specified time period
    if (time < startTime || (time-startTime) % timeInc != 0) {
      // Skip the data in irrelevant time frames
      for (int i = 0; i < npoints; i++) {
	getline(vecReader, line);
      }
    } else if (time > endTime) {
      break;
    } else {
      ibin = static_cast<int>((time-startTime)/timeInc);
      // Compute the pairwise alignment order
      w2avg = 0.0;
      for (int i = 0; i < npoints; i++) {
	getline(vecReader, line);
	ss.clear();
	ss.str(line);
	ss >> w >> vx >> vy;
	align[i][0] = w;
	align[i][1] = atan2(vy,vx);
	w2avg += w*w;
      }
      w2avg /= static_cast<double>(npoints);
      for (int i = 0; i < npoints; i++) {
	for (int j = 0; j < i; j++) {
	  dx = ddiff(lx, pos[ibin][i][0], pos[ibin][j][0]);
	  dy = ddiff(ly, pos[ibin][i][1], pos[ibin][j][1]);
	  dr = sqrt(dx*dx + dy*dy);
	  cos2dt = cos(2.0*(align[i][1]-align[j][1]));
	  wcos2dt = align[i][0]*align[j][0]*cos2dt/w2avg;
	  if (dr >= rmin && dr < rmax) {
	    rbin = static_cast<int>((dr-rmin)/rinc);
	    alignCorr[rbin][0] += cos2dt;
	    alignCorr[rbin][1] += wcos2dt;
	    alignCount[rbin]++;
	  }
	}
      }
    }
  }
  vecReader.close();

  // Normalise
  for (int i = 0; i < rnbins; i++) {
    if (alignCount[i] > 0) {
      alignCorr[i][0] /= static_cast<double>(alignCount[i]);
      alignCorr[i][1] /= static_cast<double>(alignCount[i]);
    }
  }

  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Error: cannot open the file " << outFile << endl;
    return 1;
  }
  double rlo, rmid, rhi;
  for (int i = 0; i < rnbins; i++) {
    rlo = rmin+i*rinc;
    rmid = rlo+rinc*0.5;
    rhi = rlo+rinc;
    writer << rlo << " " << rmid << " " << rhi << " "
	   << alignCorr[i][0] << " " << alignCorr[i][1] << "\n";
  }
  writer.close();

  // Clean up
  deleteArray(pos);
  deleteArray(align);
  deleteArray(alignCorr);
  deleteArray(alignCount);
}

double sgnd(double val) {
  return (0.0 < val) - (val < 0.0);
}

double ddiff(double len, double d1, double d2) {
  double dd1 = d1-d2;
  double add1 = fabs(dd1);
  double add2 = len-add1;
  return add1 < add2 ? dd1 : -sgnd(dd1)*add2;
}
