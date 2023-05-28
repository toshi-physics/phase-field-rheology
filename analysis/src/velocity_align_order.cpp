// velocity_align_order.cpp
// A program which computes the global alignment order of the system
// (both polar (Vicsek) and nematic order)

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
  if (argc != 8) {
    cout << "Usage: velocity_align_order npoints v0 startTime endTime timeInc "
	 << "velFile outFile" << endl;
    return 1;
  }

  int argi = 0;
  int npoints = stoi(string(argv[++argi]), nullptr, 10);
  double v0 = stod(string(argv[++argi]), nullptr);
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
  double vax, vay, vcmx, vcmy, v, vx, vy, w;
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
	ss >> vax >> vay >> vcmx >> vcmy;
	v = sqrt(vcmx*vcmx+vcmy*vcmy);
	vx = vcmx/v;
	vy = vcmy/v;
	w = v/v0;
	vxavg += vx;
	vyavg += vy;
	wvxavg += w*vx;
	wvyavg += w*vy;
	cos2tavg += (vx*vx-vy*vy);
	sin2tavg += (2.0*vx*vy);
	wcos2tavg += w*(vx*vx-vy*vy);
	wsin2tavg += w*(2.0*vx*vy);
	wavg += w;
      }
      vxavg /= static_cast<double>(npoints);
      vyavg /= static_cast<double>(npoints);
      wvxavg /= static_cast<double>(npoints);
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
