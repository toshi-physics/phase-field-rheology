// align_order.cpp
// A program which computes the global polar and nematic order of the system

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::string;

int main(int argc, char* argv[]) {
  if (argc != 7) {
    cout << "Usage: align_order npoints startTime endTime timeInc "
	 << "vecFile outFile" << endl;
    return 1;
  }

  int argi = 0;
  int npoints = stoi(string(argv[++argi]), nullptr, 10);
  long startTime = stol(string(argv[++argi]), nullptr, 10);
  long endTime = stol(string(argv[++argi]), nullptr, 10);
  long timeInc = stol(string(argv[++argi]), nullptr, 10);
  string vecFile (argv[++argi]);
  string outFile (argv[++argi]);

  ifstream reader;
  reader.open(vecFile);
  if (!reader) {
    cout << "Error: cannot open the file " << vecFile << endl;
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
  double vx, vy, v, w, cos2t, sin2t;
  double vxavg, vyavg, cos2tavg, sin2tavg;
  double wvxavg, wvyavg, wcos2tavg, wsin2tavg, wavg;
  double polarOrder, nematicOrder, weightedPolarOrder, weightedNematicOrder;
  double n = static_cast<double>(npoints);
  while (getline(reader, line)) {
    // Read the two header lines and get time
    getline(reader, line);
    ss.clear();
    ss.str(line);
    ss >> str >> time;

    // Only use the data from the specified time period
    if (time < startTime || (time-startTime) % timeInc != 0) {
      // Skip the data in irrelevant time frames
      for (int i = 0; i < npoints; i++) {
	getline(reader, line);
      }
    } else if (time > endTime) {
      break;
    } else {
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
	ss >> w >> vx >> vy;
	v = sqrt(vx*vx+vy*vy);
	vx /= v;
	vy /= v;
	vxavg += vx;
	vyavg += vy;
	wvxavg += w*vx;
	wvyavg += w*vy;
	cos2t = vx*vx-vy*vy;
	sin2t = 2*vx*vy;
	cos2tavg += cos2t;
	sin2tavg += sin2t;
	wcos2tavg += w*cos2t;
	wsin2tavg += w*sin2t;
	wavg += w;
      }
      // Normalise
      vxavg /= n;
      vyavg /= n;
      wvxavg /= n;
      wvyavg /= n;
      cos2tavg /= n;
      sin2tavg /= n;
      wcos2tavg /= n;
      wsin2tavg /= n;
      wavg /= n;
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
