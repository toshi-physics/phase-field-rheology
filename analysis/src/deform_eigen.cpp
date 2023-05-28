// deform_eigen.cpp
// A program which computes the eigenvalues and eigenvectors of the
// deformation (or structural) tensor

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <armadillo>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::string;
using std::vector;
using namespace arma;

struct EigenData {
  mat matrix = mat(2,2);
  vec eigenVal = vec(2);
  mat eigenVec = mat(2,2);
};

void eigen(double gxx, double gyy, double gxy, EigenData& data);

int main(int argc, char* argv[]) {
  
  if (argc != 7) {
    cout << "Usage: gyr_eigen npoints startTime endTime timeInc "
	 << "deformFile outFile" << endl;
    return 1;
  }
  
  int argi {};
  int npoints = stoi(string(argv[++argi]), nullptr, 10);
  long startTime = stol(string(argv[++argi]), nullptr, 10);
  long endTime = stol(string(argv[++argi]), nullptr, 10);
  long timeInc = stol(string(argv[++argi]), nullptr, 10);
  string deformFile (argv[++argi]);
  string outFile (argv[++argi]);
  
  ifstream reader;
  reader.open(deformFile);
  if (!reader) {
    cout << "Error: cannot open the file " << deformFile << endl;
    return 1;
  }
  
  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Error: cannot open the file " << outFile << endl;
    return 1;
  }
  
  string line, str;
  stringstream ss;
  long time; 
  double gxx, gyy, gxy;
  EigenData data;
  while (getline(reader, line)) {
    // Read the two header lines and get time
    getline(reader, line);
    ss.clear();
    ss.str(line);
    ss >> str >> time;
    
    // Only use the data from the specified time period
    if (time < startTime || (time-startTime) % timeInc != 0) {
      // Skip the data in that time frame
      for (int i = 0; i < npoints; i++) {
	getline(reader, line);
      }
    } else if (time > endTime) {
      break;
    } else {
      // Compute eigenvalues and eigenvectors
      writer << "Cells: " << npoints << "\n";
      writer << "Timestep: " << time << "\n";
      for (int i = 0; i < npoints; i++) {
	getline(reader, line);
	ss.clear();
	ss.str(line);
	ss >> gxx >> gyy >> gxy;
	eigen(gxx, gyy, gxy, data);
	// Output the eigenvalues for both trace and traceless version of
	// the tensor and the eigenvectors
	writer << data.eigenVal.at(0) << " " 
	       << 0.5*(data.eigenVal.at(0)-data.eigenVal.at(1)) << " "
	       << data.eigenVec.at(0,0) << " " << data.eigenVec.at(1,0) << " " 
	       << data.eigenVal.at(1) << " "
	       << 0.5*(data.eigenVal.at(1)-data.eigenVal.at(0)) << " "
	       << data.eigenVec.at(0,1) << " " << data.eigenVec.at(1,1) 
	       << "\n";
      }
    }
  }
  reader.close();
  writer.close();
}

void eigen(double gxx, double gyy, double gxy, EigenData& data) {
  data.matrix.at(0,0) = gxx;
  data.matrix.at(0,1) = gxy;
  data.matrix.at(1,0) = gxy;
  data.matrix.at(1,1) = gyy;
  eig_sym(data.eigenVal, data.eigenVec, data.matrix);
}
