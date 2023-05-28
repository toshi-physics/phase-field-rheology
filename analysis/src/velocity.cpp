// velocity.cpp
// A program to compute the average velocity/forces acting on the cells

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

struct Obs {
  int count = 0;
  double avg = 0.0;
  double avgSq = 0.0;
  double var = 0.0;
  double stdev = 0.0;
  double stderr = 0.0;
  void addData(double val);
  void computeStats();
};

void readObs(stringstream& ss, Obs& obs);

int main(int argc, char* argv[]) {
  if (argc != 6) {
    cout << "Usage: velocity npoints startTime endTime "
	 << "velFile outFile" << endl;
    return 1;
  }

  int argi = 0;
  int npoints = stoi(string(argv[++argi]), nullptr, 10);
  long startTime = stol(string(argv[++argi]), nullptr, 10);
  long endTime = stol(string(argv[++argi]), nullptr, 10);
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
  
  while (getline(reader, line)) {
    // Read the two header lines and get time
    getline(reader, line);
    ss.clear();
    ss.str(line);
    ss >> str >> time;
    
    // Only use the data from the specified time period
    if (time < startTime) {
      // Skip the data in that time frame
      for (int i = 0; i < npoints; i++) {
	getline(reader, line);
      }
    } else if (time > endTime) {
      break;
    } else {
      // Compute velocities
      Obs vadv, vcm, vfrict, vact;
      for (int i = 0; i < npoints; i++) {
	getline(reader, line);
	ss.clear();
	ss.str(line);
	readObs(ss, vadv);
	readObs(ss, vcm);
	readObs(ss, vfrict);
	readObs(ss, vact);
      }
      
      // Output results
      vadv.computeStats();
      vcm.computeStats();
      vfrict.computeStats();
      vact.computeStats();
      writer << time << " ";
      writer << vadv.avg << " " << vadv.stdev << " " << vadv.stderr << " ";
      writer << vcm.avg << " " << vcm.stdev << " " << vcm.stderr << " ";
      writer << vfrict.avg << " " << vfrict.stdev << " " 
	     << vfrict.stderr << " ";
      writer << vact.avg << " " << vact.stdev << " " << vact.stderr << "\n";
    }
  }
  reader.close();
  writer.close();
}

void readObs(stringstream& ss, Obs& obs) {
  double vx, vy, v;
  ss >> vx >> vy;
  v = sqrt(vx*vx+vy*vy);
  obs.addData(v);
}

// Member functions of Obs
void Obs::addData(double val) {
  avg += val;
  avgSq += val*val;
  count++;
}

void Obs::computeStats() {
  double n = static_cast<double>(count);
  if (count <= 1) return;
  avg /= n;
  avgSq /= n;
  var = n/(n-1.0)*(avgSq-avg*avg);
  stdev = sqrt(var);
  stderr = stdev/sqrt(n);
}
