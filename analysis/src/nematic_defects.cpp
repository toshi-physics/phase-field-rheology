// nematic_defects.cpp
// A program which identifies nematic defects of the Q tensor field

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <armadillo>
#include "array.hpp"

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::string;
using std::vector;
//using namespace arma;

// Helper functions
void smooth(int lx, int ly, int binsize, double*** field, double*** smoothed);
vector<vector<int> > findMinima(int lx, int ly, int dist, double** field);
double defectCharge(int lx, int ly, int i, int j, double*** dirField,
		    const vector<vector<int> >& contour);
bool printField(int lx, int ly, double** field, string file);
bool printVecField(int lx, int ly, double*** field, string file);
int iup(int len, int i);
int idown(int len, int i);
int iwrap(int len, int i);
double cgrad4(int i, int j, int uu, int u, int d, int dd, int ic, int oc,
	      double*** field);

struct Defect {
  double charge;
  int x;
  int y;
};

int main(int argc, char* argv[]) {
  if (argc != 11) {
    cout << "Usage: nematic_defects lx ly binsize minDist contourDist "
	 << "startTime endTime timeInc fieldFileName outFile" << endl;
    return 1;
  }

  int argi = 0;
  int lx = stoi(string(argv[++argi]), nullptr, 10);
  int ly = stoi(string(argv[++argi]), nullptr, 10);
  int binsize = stoi(string(argv[++argi]), nullptr, 10);
  int minDist = stoi(string(argv[++argi]), nullptr, 10);
  int contourDist = stoi(string(argv[++argi]), nullptr, 10);
  long startTime = stol(string(argv[++argi]), nullptr, 10);
  long endTime = stol(string(argv[++argi]), nullptr, 10);
  long timeInc = stol(string(argv[++argi]), nullptr, 10);
  string fieldFile (argv[++argi]);
  string outFile (argv[++argi]);

  ifstream reader;
  string fullFieldFile, line;
  stringstream ss;
  int x, y;
  double qxx, qxy;
  double pi = M_PI;
  double*** field = create3DArray<double>(lx, ly, 2);
  double*** smoothedField = create3DArray<double>(lx, ly, 2);
  double** topoCharge = create2DArray<double>(lx, ly);
  double** eigenVal = create2DArray<double>(lx, ly);
  double*** eigenVec = create3DArray<double>(lx, ly, 2);
  const double minDefectCharge = 0.495;

  // A closed square contour for summing angles
  vector<vector<int> > contour;
  for (int i = -contourDist; i < contourDist; i++) {
    contour.push_back({contourDist,i});
  }
  for (int i = contourDist; i > -contourDist; i--) {
    contour.push_back({i,contourDist});
  }
  for (int i = contourDist; i > -contourDist; i--) {
    contour.push_back({-contourDist,i});
  }
  for (int i = -contourDist; i < contourDist; i++) {
    contour.push_back({i,-contourDist});
  }
  
  for (long time = startTime; time <= endTime; time += timeInc) {
    ss.clear();
    ss << fieldFile << "." << time;
    fullFieldFile = ss.str();
    reader.open(fullFieldFile);
    if (!reader) {
      cout << "Error: cannot open the file " << fullFieldFile << endl;
      return 1;
    }
    while (getline(reader, line)) {
      ss.clear();
      ss.str(line);
      ss >> x >> y >> qxx >> qxy;
      field[x][y][0] = qxx; // Qxx = -Qyy
      field[x][y][1] = qxy;
    }
    
    // Smooth the field by a square window
    smooth(lx, ly, binsize, field, smoothedField);
    
    // Computing topological charge density
    int iu, iuu, id, idd, ju, juu, jd, jdd;
    for (int i = 0; i < lx; i++) {
      iu = iup(lx,i);
      iuu = iup(lx,iu);
      id = idown(lx,i);
      idd = idown(lx,id);
      for (int j = 0; j < ly; j++) {
	ju = iup(ly,j);
	juu = iup(ly,ju);
	jd = idown(ly,j);
	jdd = idown(ly,jd);
	topoCharge[i][j] = 1.0/(2.0*pi) * 
	  (cgrad4(i,j,iuu,iu,id,idd,0,0,smoothedField) * 
	   cgrad4(i,j,juu,ju,jd,jdd,1,1,smoothedField) -
	   cgrad4(i,j,iuu,iu,id,idd,1,0,smoothedField) *
	   cgrad4(i,j,juu,ju,jd,jdd,0,1,smoothedField));
      }
    }
    printField(lx, ly, topoCharge, "topo_charge.dat");

    // Compute eigenvalues
    for (int i = 0; i < lx; i++) {
      for (int j = 0; j < ly; j++) {
	eigenVal[i][j] = sqrt(smoothedField[i][j][0] * smoothedField[i][j][0] +
			      smoothedField[i][j][1] * smoothedField[i][j][1]);
      }
    }
    printField(lx, ly, eigenVal, "eigenvalue.dat");

    // Compute major deform axis (eigenvector with +s)
    double v, vx, vy;
    /*mat Q (2,2);
    mat evec (2,2);
    vec eval (2);*/
    for (int i = 0; i < lx; i++) {
      for (int j = 0; j < ly; j++) {
	vx = 1.0;
	vy = (eigenVal[i][j] - smoothedField[i][j][0]) / 
	  smoothedField[i][j][1];
	v = sqrt(vx*vx+vy*vy);
	eigenVec[i][j][0] = vx / v;
	eigenVec[i][j][1] = vy / v;
	/*Q(0,0) = smoothedField[i][j][0];
	Q(0,1) = smoothedField[i][j][1];
	Q(1,0) = smoothedField[i][j][1];
	Q(1,1) = -smoothedField[i][j][0];
	eig_sym(eval, evec, Q);
	eigenVec[i][j][0] = evec(1,0);
	eigenVec[i][j][1] = evec(1,1);*/
      }
    }
    printVecField(lx, ly, eigenVec, "eigenvector.dat");

    // Find defects (find minima in eigenvalue field)
    vector<vector<int> > defectpts = findMinima(lx, ly, minDist, eigenVal);
    int npts = static_cast<int>(defectpts.size());
    ofstream writer;
    writer.open("defects.dat");
    double q;
    for (int i = 0; i < npts; i++) {
      x = defectpts[i][0];
      y = defectpts[i][1];
      q = defectCharge(lx, ly, x, y, eigenVec, contour);
      if (q < minDefectCharge) continue;
      writer << x << " " << y << " "
	     << (topoCharge[x][y] > 0.0 ? 1 : -1) << "\n";
    }
    writer.close();
  }

  // Clean up
  deleteArray(field);
  deleteArray(smoothedField);
  deleteArray(topoCharge);
  deleteArray(eigenVal);
  deleteArray(eigenVec);
}

void smooth(int lx, int ly, int binsize, double*** field, double*** smoothed) {
  int start = -binsize/2;
  int end = start+binsize;
  int n = binsize*binsize;
  cout << start << " " << end << endl;
#pragma omp parallel for default(none)		\
  shared(lx, ly, start, end, n, field, smoothed)
  for (int i = 0; i < lx; i++) {
    for (int j = 0; j < ly; j++) {
      double sumxx = 0.0; // Qxx = -Qyy
      double sumxy = 0.0;
      int ki, lj;
      for (int k = start; k < end; k++) {
	ki = iwrap(lx, i+k);
	for (int l = start; l < end; l++) {
	  lj = iwrap(ly, j+l);
	  sumxx += field[ki][lj][0];
	  sumxy += field[ki][lj][1];
	}
      }
      smoothed[i][j][0] = sumxx / n;
      smoothed[i][j][1] = sumxy / n;
    }
  }
}

vector<vector<int> > findMinima(int lx, int ly, int dist, double** field) {
  int binsize = dist*2+1;
  vector<int> pt (2);
  vector<vector<int> > pts;
  int start = -binsize/2;
  int end = binsize+start;
  bool smaller;
  for (int i = 0; i < lx; i++) {
    for (int j = 0; j < ly; j++) {
      smaller = true;
      int ki, lj;
      for (int k = start; k < end && smaller; k++) {
	ki = iwrap(lx,i+k);
	for (int l = start; l < end && smaller; l++) {
	  lj = iwrap(ly,j+l);
	  if (ki == i && lj == j) continue;
	  smaller = field[i][j] < field[ki][lj];
	}
      }
      if (smaller) {
	pt[0] = i;
	pt[1] = j;
	pts.push_back(pt);
      }
    }
  }
  return pts;
}

double defectCharge(int lx, int ly, int i, int j, double*** dirField,
		    const vector<vector<int> >& contour) {
  double pi = M_PI;
  double twopi = pi*2.0;
  double halfpi = pi/2.0;
  double angleSum = 0.0;
  double angle;
  int ic, jc, ip, jp, kup;
  int n = contour.size();
  for (int k = 0; k < n; k++) {
    kup = iup(n,k);
    ic = iwrap(lx, i+contour[kup][0]);
    jc = iwrap(ly, j+contour[kup][1]);
    ip = iwrap(lx, i+contour[k][0]);
    jp = iwrap(ly, j+contour[k][1]);
    angle = acos(dirField[ic][jc][0]*dirField[ip][jp][0]+
		 dirField[ic][jc][1]*dirField[ip][jp][1]);
    //cout << i << " " << j << " " << ic << " " << jc << " " << angle << endl;
    //angle = angle > halfpi ? pi - angle : angle;
    //cout << i << " " << j << " " << ic << " " << jc << " " << angle << endl;
    //angleSum += angle;
    angleSum += angle > halfpi ? pi - angle : angle;
  }
  return angleSum / twopi;
}

bool printField(int lx, int ly, double** field, string file) {
  ofstream writer;
  writer.open(file);
  if (!writer) {
    cout << "Error: cannot open the file " << file << endl;
    return false;
  }
  for (int i = 0; i < lx; i++) {
    for (int j = 0; j < ly; j++) {
      writer << i << " " << j << " " << field[i][j] << "\n";
    }
    writer << "\n";
  }
  writer.close();
  return true;
}

bool printVecField(int lx, int ly, double*** field, string file) {
  ofstream writer;
  writer.open(file);
  if (!writer) {
    cout << "Error: cannot open the file " << file << endl;
    return false;
  }
  for (int i = 0; i < lx; i++) {
    for (int j = 0; j < ly; j++) {
      writer << i << " " << j << " " 
	     << field[i][j][0] << " " << field[i][j][1] << "\n";
    }
    writer << "\n";
  }
  writer.close();
  return true;
}


inline int iup(int len, int i) {
  return (i+1 >= len) ? 0 : i+1;
}

inline int idown(int len, int i) {
  return (i-1 < 0) ? len-1 : i-1;
}

inline int iwrap(int len, int i) {
  int remainder = i % len;
  return remainder >= 0 ? remainder : len + remainder;
}

inline double cgrad4(int i, int j, int uu, int u, int d, int dd, int ic,
                     int oc, double*** field) {
  switch (oc) {
  case 0: return (-field[uu][j][ic] +
                  8.0 * (field[u][j][ic] - field[d][j][ic]) +
                  field[dd][j][ic]) / 12.0;
  case 1: return (-field[i][uu][ic] +
                  8.0 * (field[i][u][ic] - field[i][d][ic]) +
                  field[i][dd][ic]) / 12.0;
  default: return 0.0;
  }
}
