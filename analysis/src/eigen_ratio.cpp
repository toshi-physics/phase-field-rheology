// eigen_align_order.cpp
// A program which computes the global alignment order of the system
// (both polar (Vicsek) and nematic order) based on given eigenvalues and
// eigenvectors

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::string;
using std::vector;

int main(int argc, char* argv[]) {
  if (argc != 7) {
    cout << "Usage: eigen_ratio npoints startTime endTime timeInc "
	 << "eigenFile outFile" << endl;
    return 1;
  }

  int argi = 0;
  int npoints = stoi(string(argv[++argi]), nullptr, 10);
  long startTime = stol(string(argv[++argi]), nullptr, 10);
  long endTime = stol(string(argv[++argi]), nullptr, 10);
  long timeInc = stol(string(argv[++argi]), nullptr, 10);
  string eigenFile (argv[++argi]);
  string outFile (argv[++argi]);

  ifstream reader;
  reader.open(eigenFile);
  if (!reader) {
    cout << "Error: cannot open the file " << eigenFile << endl;
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
  double lam1, tlam1, vx1, vy1, lam2, tlam2, vx2, vy2, ratio;
  double ratioavg;
  
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
      ratioavg = 0.0;
      for (int i = 0; i < npoints; i++) {
	getline(reader, line);
	ss.clear();
	ss.str(line);
	// tlam1 and tlam2 are eigenvalues corresponding to the traceless
	// version of the tensor
	// Assuming the eigevalues are sorted in ascending order
	// (i.e., the second one is the larger eigenvalue)
	ss >> lam1 >> tlam1 >> vx1 >> vy1 >> lam2 >> tlam2 >> vx2 >> vy2;
	ratio = fabs(lam2-lam1)/(lam1+lam2);
	ratioavg += ratio;
      }
      ratioavg /= static_cast<double>(npoints);
      writer << time << " " << ratioavg << " \n";
    }
  }
  reader.close();
  writer.close();
}
