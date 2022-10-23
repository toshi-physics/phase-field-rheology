// vicsek_order.cpp
// A program which computes the Vicsek order parameter of the system

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
    cout << "Usage: velocity npoints startTime endTime timeInc "
	 << "velFile outFile" << endl;
    return 1;
  }

  int argi = 0;
  int npoints = stoi(string(argv[++argi]), nullptr, 10);
  long startTime = stol(string(argv[++argi]), nullptr, 10);
  long endTime = stol(string(argv[++argi]), nullptr, 10);
  long timeInc = stol(string(argv[++argi]), nullptr, 10);
  string velFile (argv[++argi]);
  string outFile (argv[++argi]);

  ifstream reader;
  reader.open(velFile);
  if (!reader) {
    cout << "Error: cannot open the file " << velFile << endl;
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
  double vax, vay, vcmx, vcmy, v, vxavg, vyavg;
  double order;
  
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
      for (int i = 0; i < npoints; i++) {
	getline(reader, line);
	ss.clear();
	ss.str(line);
	ss >> vax >> vay >> vcmx >> vcmy;
	v = sqrt(vcmx*vcmx+vcmy*vcmy);
	vxavg += vcmx/v;
	vyavg += vcmy/v;
      }
      vxavg /= static_cast<double>(npoints);
      vyavg /= static_cast<double>(npoints);
      order = sqrt(vxavg*vxavg+vyavg*vyavg);
      writer << time << " " << order << "\n";
    }
  }
  reader.close();
  writer.close();
}
