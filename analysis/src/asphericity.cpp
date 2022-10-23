// asphericity.cpp
// A program to compute the average asphericity

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::string;
using std::vector;

int main(int argc, char* argv[]) {
  
  if (argc != 6) {
    cout << "usage: asphericity npoints startTime endTime "
	 << "shapeFile outFile" << endl;
    return 1;
  }

  int argi {};
  int npoints {stoi(string(argv[++argi]), nullptr, 10)};
  long startTime {stoi(string(argv[++argi]), nullptr, 10)};
  long endTime {stoi(string(argv[++argi]), nullptr, 10)};
  string shapeFile {argv[++argi]};
  string outFile {argv[++argi]};
  
  ifstream reader;
  reader.open(shapeFile);
  if (!reader) {
    cout << "Problem with opening gyration file!" << endl;
    return 1;
  }

  ofstream writer;
  writer.open(outFile);
  if (!writer) {
    cout << "Problem with opening output file!" << endl;
    return 1;
  }
  writer << std::setprecision(10) << std::fixed;
  
  string line, str;
  istringstream iss;
  long time; 
  int pixels;
  double area, perimeter, pixelArea, chainPerimeter, asphere;
  double asphereAvg {};
  double asphereAvgSq {};
  double n = static_cast<double>(npoints);
  while (getline(reader, line)) {
    // Read the two header lines and get time
    getline(reader, line);
    iss.clear();
    iss.str(line);
    iss >> str >> time;
    
    // Only use the data from the specified time period
    if (time < startTime) {
      // Skip the data in that time frame
      for (int i {}; i < npoints; i++) {
	getline(reader, line);
      }
    } else if (time > endTime) {
      break;
    } else {
      // Compute asphericity
      asphereAvg = 0.0;
      asphereAvgSq = 0.0;
      for (int i {}; i < npoints; i++) {
	getline(reader, line);
	iss.clear();
	iss.str(line);
	iss >> perimeter >> area >> chainPerimeter >> pixelArea >> pixels;
	asphere = perimeter*perimeter/(4.0*M_PI*area);
	asphereAvg += asphere;
	asphereAvgSq += asphere*asphere;
      }

      // Normalise
      asphereAvg /= n;
      asphereAvgSq /= n;
      double var {n/(n-1.0)*(asphereAvgSq-asphereAvg*asphereAvg)};
      double stdev {sqrt(var)};
      double stderr {stdev/sqrt(n)};
      
      // Output results to file
      writer << time << " " << std::setprecision(10) << std::fixed 
	     << asphereAvg << " " << stdev << " " << stderr << endl;
      writer.unsetf(std::ios_base::floatfield);
    }
  }
  reader.close();
  writer.close();
}
