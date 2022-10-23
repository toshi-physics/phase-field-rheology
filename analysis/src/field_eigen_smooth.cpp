// field_eigen.cpp
// A code which computes the eigenvalues and eigenvectors at each point
// of a tensor field

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <armadillo>
#include "array.hpp"

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::string;
using namespace arma;

struct EigenData {
  mat matrix = mat(2,2);
  vec eigenVal = vec(2);
  mat eigenVec = mat(2,2);
};

void eigen(double gxx, double gyy, double gxy, EigenData& data);
double len(double* x);
int wrap(int lx, int i);

int main(int argc, char* argv[]) {
  if (argc != 7) {
    cout << "Usage: field_eigen_smooth lx ly width smoothType " 
	 << "fieldFile outFile" << endl;
    return 1;
  }

  int argi = 0;
  int lx = stoi(string(argv[++argi]), nullptr, 10);
  int ly = stoi(string(argv[++argi]), nullptr, 10);
  int width = stoi(string(argv[++argi]), nullptr, 10);
  string smoothType (argv[++argi]);
  string fieldFile (argv[++argi]);
  string outFile (argv[++argi]);

  double*** field = create3DArray<double>(lx, ly, 3);
  double*** smoothField = create3DArray<double>(lx, ly, 2);

  ifstream reader;
  reader.open(fieldFile);
  if (!reader) {
    cout << "Error: cannot open the file " << fieldFile << endl;
    return 1;
  }

  string line;
  stringstream ss;
  int x, y;
  double txx, txy, tyy, e0, e1, ratio;
  EigenData data;
  while (getline(reader, line)) {
    if (line.size() == 0) continue;
    ss.clear();
    ss.str(line);
    ss >> x >> y >> txx >> tyy >> txy;
    eigen(txx, tyy, txy, data);
    e0 = data.eigenVal.at(0);
    e1 = data.eigenVal.at(1);
    ratio = (e1-e0)/(e1+e0);
    field[x][y][0] = ratio*data.eigenVec.at(0,1);
    field[x][y][1] = ratio*data.eigenVec.at(1,1);
    field[x][y][2] = ratio;
  }

  // Do smoothing
  for (int i = 0; i < lx; i++) {
    for (int j = 0; j < ly; j++) {
      for (int k = 0; k < width; k++) {
	x = wrap(lx,i-width/2+k);
	for (int l = 0; l < width; l++) {
	  y = wrap(ly,j-width/2+l);
	  for (int m = 0; m < 2; m++) {
	    smoothField[i][j][m] += field[x][y][m];
	  }
	}
      }
      for (int m = 0; m < 2; m++) {
	smoothField[i][j][m] /= static_cast<double>(width*width);
      }
    }
  }

  // Output results
  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Error: cannot open the file " << outFile << endl;
    return 1;
  }

  double vx, vy, r;
  for (int i = 0; i < lx; i++) {
    for (int j = 0; j < ly; j++) {
      r = len(smoothField[i][j]);
      vx = smoothField[i][j][0]/r;
      vy = smoothField[i][j][1]/r;
      writer << i << " " << j << " " << field[i][j][0] << " " 
	     << field[i][j][1] << " " << field[i][j][2] << " "
	     << vx << " " << vy << " " << r << "\n";
    }
  }

  // Clean up
  deleteArray(field);
  deleteArray(smoothField);
}

void eigen(double gxx, double gyy, double gxy, EigenData& data) {
  data.matrix.at(0,0) = gxx;
  data.matrix.at(0,1) = gxy;
  data.matrix.at(1,0) = gxy;
  data.matrix.at(1,1) = gyy;
  eig_sym(data.eigenVal, data.eigenVec, data.matrix);
}

inline double len(double* x) {
  double sum = 0.0;
  for (int i = 0; i < 2; i++) {
    sum += x[i]*x[i];
  }
  return sqrt(sum);
}

inline int wrap(int lx, int i) {
  int remainder = i % lx;
  if (remainder >= 0) return remainder;
  return lx + remainder;
}
