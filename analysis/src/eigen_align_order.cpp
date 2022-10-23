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
    cout << "Usage: eigen_align_order npoints startTime endTime timeInc "
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
  double lam1, tlam1, vx1, vy1, lam2, tlam2, vx2, vy2, w;
  double vxavg, vyavg, cos2tavg, sin2tavg, wavg;
  double wvxavg, wvyavg, wcos2tavg, wsin2tavg;
  double polarOrder, nematicOrder, weightedPolarOrder, weightedNematicOrder;
  
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
      // Compute velocities
      vxavg = 0.0;
      vyavg = 0.0;
      wvxavg = 0.0;
      wvyavg = 0.0;
      cos2tavg = 0.0;
      sin2tavg = 0.0;
      wcos2tavg = 0.0;
      wsin2tavg = 0.0;
      wavg = 0.0;
      for (int i = 0; i < npoints; i++) {
	getline(reader, line);
	ss.clear();
	ss.str(line);
	// tlam1 and tlam2 are eigenvalues corresponding to the traceless
	// version of the tensor
	// Assuming the eigevalues are sorted in ascending order
	// (i.e., the second one is the larger eigenvalue)
	ss >> lam1 >> tlam1 >> vx1 >> vy1 >> lam2 >> tlam2 >> vx2 >> vy2;
	w = fabs(lam2-lam1)/(lam1+lam2);
	vxavg += vx2;
	vyavg += vy2;
	wvxavg += w*vx2;
	wvyavg += w*vy2;
	cos2tavg += (vx2*vx2-vy2*vy2);
	sin2tavg += (2.0*vx2*vy2);
	wcos2tavg += w*(vx2*vx2-vy2*vy2);
	wsin2tavg += w*(2.0*vx2*vy2);
	wavg += w;
      }
      vxavg /= static_cast<double>(npoints);
      vyavg /= static_cast<double>(npoints);
      wvxavg /= static_cast<double>(npoints);
      wvyavg /= static_cast<double>(npoints);
      wvyavg /= static_cast<double>(npoints);
      cos2tavg /= static_cast<double>(npoints);
      sin2tavg /= static_cast<double>(npoints);
      wcos2tavg /= static_cast<double>(npoints);
      wsin2tavg /= static_cast<double>(npoints);
      wavg /= static_cast<double>(npoints);
      polarOrder = sqrt(vxavg*vxavg+vyavg*vyavg);
      nematicOrder = sqrt(cos2tavg*cos2tavg+sin2tavg*sin2tavg);
      weightedPolarOrder = sqrt(wvxavg*wvxavg+wvyavg*wvyavg);
      weightedNematicOrder = sqrt(wcos2tavg*wcos2tavg+wsin2tavg*wsin2tavg);
      writer << time << " " << wavg << " " << polarOrder << " " 
	     << weightedPolarOrder << " " << weightedPolarOrder/wavg << " "
	     << nematicOrder << " " << weightedNematicOrder << " "
	     << weightedNematicOrder/wavg << "\n";
    }
  }
  reader.close();
  writer.close();
}
