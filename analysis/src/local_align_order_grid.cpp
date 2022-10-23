// local_align_order.cpp
// A program which computes the local nematic order of the system

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

double distSq(double* len, double* x1, double* x2);
double sgnd(double val);
double ddiff(double len, double d1, double d2);
void error(int n, double avg, double avgSq, double& dev, double& err);

int main(int argc, char* argv[]) {
  if (argc != 15) {
    cout << "Usage: local_align_order npoints lx ly ngridx ngridy "
	 << "minRadius maxRadius radiusInc startTime endTime timeInc "
	 << "posFile vecFile outFile" << endl;
    return 1;
  }

  double boxsize[2];
  int argi = 0;
  int npoints = stoi(string(argv[++argi]), nullptr, 10);
  boxsize[0] = stod(string(argv[++argi]), nullptr);
  boxsize[1] = stod(string(argv[++argi]), nullptr);
  int ngridx = stoi(string(argv[++argi]), nullptr, 10);
  int ngridy = stoi(string(argv[++argi]), nullptr, 10);
  double minRadius = stod(string(argv[++argi]), nullptr);
  double maxRadius = stod(string(argv[++argi]), nullptr);
  double radiusInc = stod(string(argv[++argi]), nullptr);
  long startTime = stol(string(argv[++argi]), nullptr, 10);
  long endTime = stol(string(argv[++argi]), nullptr, 10);
  long timeInc = stol(string(argv[++argi]), nullptr, 10);
  string posFile (argv[++argi]);
  string vecFile (argv[++argi]);
  string outFile (argv[++argi]);

  // Read position data
  PositionReader posReader;
  if (!posReader.open(posFile, npoints, boxsize[0], boxsize[1], timeInc)) {
    cout << "Error: cannot open the file " << posFile << endl;
    return 1;
  }

  ifstream vecReader;
  vecReader.open(vecFile);
  if (!vecReader) {
    cout << "Error: cannot open the file " << vecFile << endl;
    return 1;
  }
  
  bool foundPosData, foundAlignData;
  const int navg = 9;
  const int norders = 4;
  int ntbins = static_cast<int>((endTime-startTime)/timeInc)+1;
  int nrbins = static_cast<int>((maxRadius-minRadius)/radiusInc)+1;

  double** vecData = create2DArray<double>(npoints, navg);
  double*** vecAvg = create3DArray<double>(ngridx, ngridy, navg);
  double** orderAvg = create2DArray<double>(nrbins, norders);
  double** orderAvgSq = create2DArray<double>(nrbins, norders);
  double*** distSqMat = create3DArray<double>(ngridx, ngridy, npoints);
  double* gridx = create1DArray<double>(ngridx);
  double* gridy = create1DArray<double>(ngridy);
  int** count = create2DArray<int>(ngridx, ngridy);
  double radius, radiusSq, drsq, w, v, vx, vy, wvx, wvy;
  double cos2t, sin2t, wcos2t, wsin2t;
  double polarOrder, nematicOrder, weightedPolarOrder, weightedNematicOrder;
  long t;
  stringstream ss;
  string line, str;
  
  // Compute grid positions
  double gridxWidth = boxsize[0] / static_cast<double>(ngridx);
  double gridyWidth = boxsize[1] / static_cast<double>(ngridy);
  for (int i = 0; i < ngridx; i++) {
    gridx[i] = (i+0.5)*gridxWidth;
  }
  for (int i = 0; i < ngridy; i++) {
    gridy[i] = (i+0.5)*gridyWidth;
  }

  for (long time = startTime; time <= endTime; time += timeInc) {
    // Get position data for this time frame
    foundPosData = false;
    foundAlignData = false;
    while (posReader.nextFrame()) {
      t = posReader.getTime();
      if (t == time) {
	foundPosData = true;
	break;
      } else if (t > time) {
	foundPosData = false;
	break;
      }
    }
    if (!foundPosData) {
      cout << "Error: cannot find the position data for time = " 
	   << time << endl;
      return 1;
    }
    
    while (getline(vecReader, line)) {
      // Read the two header lines and get time
      getline(vecReader, line);
      ss.clear();
      ss.str(line);
      ss >> str >> t;
      if (t < time) {
	// Skip the data in irrelevant time frames
	for (int i = 0; i < npoints; i++) {
	  getline(vecReader, line);
	}	
      } else if (t == time) {
	for (int i = 0; i < npoints; i++) {
	  getline(vecReader, line);
	  ss.clear();
	  ss.str(line);
	  ss >> w >> vx >> vy;
	  v = sqrt(vx*vx+vy*vy);
	  vx /= v;
	  vy /= v;
	  vecData[i][0] = w; // weight
	  vecData[i][1] = vx; // cos(t)
	  vecData[i][2] = vy; // sin(t)
	  vecData[i][3] = w*vx; // w*cos(t)
	  vecData[i][4] = w*vy; // w*sin(t)	  
	  vecData[i][5] = vx*vx-vy*vy; // cos(2t)
	  vecData[i][6] = 2*vx*vy; // sin(2t)
	  vecData[i][7] = w*vecData[i][5]; // w*cos(2t)
	  vecData[i][8] = w*vecData[i][6]; // w*sin(2t)
	}
	foundAlignData = true;
	break;
      } else if (t > time) {
	foundAlignData = false;
	break;
      }      
    }    
    
    if (!foundAlignData) {
      cout << "Error: cannot find the alignment data for time = " 
	   << time << endl;
      return 1;
    }

    // Compute distance matrix
    double xg[2], xp[2];
    for (int i = 0; i < ngridx; i++) {
      xg[0] = gridx[i];
      for (int j = 0; j < ngridy; j++) {
	xg[1] = gridy[j];
	for (int k = 0; k < npoints; k++) {
	  xp[0] = posReader.getPosition(k,0);
	  xp[1] = posReader.getPosition(k,1);
	  drsq = distSq(boxsize,xg,xp);
	  distSqMat[i][j][k] = drsq;
	}
      }
    }
    
    // Compute local alignment order for each radius
    //int count;
    for (int n = 0; n < nrbins; n++) {
      radius = n*radiusInc+minRadius;
      radiusSq = radius*radius;
      for (int i = 0; i < ngridx; i++) {
	for (int j = 0; j < ngridy; j++) {
	  // Reset
	  count[i][j] = 0;
	  for (int l = 0; l < navg; l++) {
	    vecAvg[i][j][l] = 0.0; 
	  }
	  for (int k = 0; k < npoints; k++) {
	    if (distSqMat[i][j][k] > radiusSq) continue;
	    for (int l = 0; l < navg; l++) {
	      vecAvg[i][j][l] += vecData[k][l];
	    }
	    count[i][j]++;
	  }
	  // Average over the cells
	  if (count[i][j] > 0) {
	    for (int l = 0; l < navg; l++) {
	      vecAvg[i][j][l] /= static_cast<double>(count[i][j]);
	    }
	  }
	}
      }
      polarOrder = 0.0;
      weightedPolarOrder = 0.0;
      nematicOrder = 0.0;
      weightedNematicOrder = 0.0;
      int ngridptsWithData = 0;
      for (int i = 0; i < ngridx; i++) {
	for (int j = 0; j < ngridy; j++) {
	  if (count[i][j] == 0) continue;
	  ngridptsWithData++;
	  w = vecAvg[i][j][0];
	  vx = vecAvg[i][j][1];
	  vy = vecAvg[i][j][2];
	  wvx = vecAvg[i][j][3];
	  wvy = vecAvg[i][j][4];
	  cos2t = vecAvg[i][j][5];
	  sin2t = vecAvg[i][j][6];
	  wcos2t = vecAvg[i][j][7];
	  wsin2t = vecAvg[i][j][8];
	  polarOrder += sqrt(vx*vx+vy*vy);
	  weightedPolarOrder += sqrt(wvx*wvx+wvy*wvy)/w;
	  nematicOrder += sqrt(cos2t*cos2t+sin2t*sin2t);
	  weightedNematicOrder += sqrt(wcos2t*wcos2t+wsin2t*wsin2t)/w;
	}
      }
      if (ngridptsWithData > 0) {
	polarOrder /= static_cast<double>(ngridptsWithData);
	weightedPolarOrder /= static_cast<double>(ngridptsWithData);
	nematicOrder /= static_cast<double>(ngridptsWithData);
	weightedNematicOrder /= static_cast<double>(ngridptsWithData);
      }
      orderAvg[n][0] += polarOrder;
      orderAvg[n][1] += weightedPolarOrder;
      orderAvg[n][2] += nematicOrder;
      orderAvg[n][3] += weightedNematicOrder;
      orderAvgSq[n][0] += polarOrder*polarOrder;
      orderAvgSq[n][1] += weightedPolarOrder*weightedPolarOrder;
      orderAvgSq[n][2] += nematicOrder*nematicOrder;
      orderAvgSq[n][3] += weightedNematicOrder*weightedNematicOrder;
    } // Close loop over radius
  } // Close loop over time

  posReader.close();
  vecReader.close();

  // Averge over time
  for (int n = 0; n < nrbins; n++) {
    for (int k = 0; k < norders; k++) {
      orderAvg[n][k] /= static_cast<double>(ntbins);
      orderAvgSq[n][k] /= static_cast<double>(ntbins);
    }
  }
  
  // Output results
  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Error: cannot open the file " << outFile << endl;
    return 1;
  }

  double stdev, stderr;
  for (int n = 0; n < nrbins; n++) {
    radius = n*radiusInc+minRadius;
    writer << radius << " ";
    for (int k = 0; k < norders; k++) {
      error(ntbins, orderAvg[n][k], orderAvgSq[n][k], stdev, stderr);
      writer << orderAvg[n][k] << " " << stdev << " " << stderr << " ";
    }
    writer << "\n";
  }
  writer.close();

  // Clean up
  deleteArray(distSqMat);
  deleteArray(gridx);
  deleteArray(gridy);
  deleteArray(vecData);
  deleteArray(vecAvg);
  deleteArray(orderAvg);
  deleteArray(orderAvgSq);
}

double distSq(double* len, double* x1, double* x2) {
  double sum = 0.0;
  double dx;
  for (int i = 0; i < 2; i++) {
    dx = ddiff(len[i],x1[i],x2[i]);
    sum += dx*dx;
  }
  return sum;
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

void error(int n, double avg, double avgSq, double& dev, double& err) {
  if (n > 1) {
    dev = sqrt(n/(n-1.0)*(avgSq-avg*avg));
    err = dev/sqrt(n);
  } else {
    dev = 0.0;
    err = 0.0;
  }
}
