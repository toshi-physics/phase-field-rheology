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

int main(int argc, char* argv[]) {
  if (argc != 3) {
    cout << "Usage: field_eigen fieldFile outFile" << endl;
    return 1;
  }

  int argi = 0;
  string fieldFile (argv[++argi]);
  string outFile (argv[++argi]);

  ifstream reader;
  reader.open(fieldFile);
  if (!reader) {
    cout << "Error: cannot open the file " << fieldFile << endl;
    return 1;
  }

  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Error: cannot open the file " << outFile << endl;
    return 1;
  }

  string line;
  stringstream ss;
  int x, y;
  double txx, txy, tyy;
  EigenData data;
  while (getline(reader, line)) {
    if (line.size() == 0) {
      writer << "\n";
      continue;
    }
    ss.clear();
    ss.str(line);
    ss >> x >> y >> txx >> tyy >> txy;
    cout << x << " " << y << " " << txx << " " << tyy << " " << txy << endl;
    eigen(txx, tyy, txy, data);
    writer << x << " " << y << " " << data.eigenVal.at(0) << " " 
	   << 0.5*(data.eigenVal.at(0)-data.eigenVal.at(1)) << " "
	   << data.eigenVec.at(0,0) << " " << data.eigenVec.at(1,0) << " " 
	   << data.eigenVal.at(1) << " "
	   << 0.5*(data.eigenVal.at(1)-data.eigenVal.at(0)) << " "
	   << data.eigenVec.at(0,1) << " " << data.eigenVec.at(1,1) 
	   << "\n";
  }
}

void eigen(double gxx, double gyy, double gxy, EigenData& data) {
  data.matrix.at(0,0) = gxx;
  data.matrix.at(0,1) = gxy;
  data.matrix.at(1,0) = gxy;
  data.matrix.at(1,1) = gyy;
  eig_sym(data.eigenVal, data.eigenVec, data.matrix);
}
