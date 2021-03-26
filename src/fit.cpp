#include <iostream>
#include <iomanip>
#include <string>

#include <TFile.h>
#include <TTree.h>

#include <cerenkovFitter.h>
#include <cerenkovHistogramFile.h>

void parameterCalculator(double * functionParameters, const double *& p, const double &age, const double &delta, const double &energy)
{
  double & t1  = functionParameters[0];
  double & t2  = functionParameters[1];
  double & nu  = functionParameters[2];
  double & eps = functionParameters[3];
  
  nu  = p[0] * std::pow(delta, p[1]) + p[2] * std::log(age);
  t1  = p[3] * std::pow(delta, p[4]) * std::pow(energy,p[11]) + p[5] * std::log(age);
  t2  = t1 * (p[6] + p[7]*(age-1));
  eps = p[8] + p[9]*std::pow(energy,p[10]);
}

int main(int argc, char ** argv)
{
  double minAge = 0.8;
  double maxAge = 1.2;
  
  std::vector<double> parameterSeed;
  parameterSeed = {0.36241, -0.101723, 1.4695, 1.23346, 0.311041, -0.0468679, 1.30104, 0.275443, 0.00340445}; // gamma

  CerenkovFitter fitter;
  fitter.SetParameterFunction(parameterCalculator);
  fitter.SetAgeInterval(minAge,maxAge);
  fitter.SetVerbosity(true);
  
  // Loop over input files
  for (int iFile = 1; iFile < argc; iFile++) {

    int maxShowers = 120;
    double weight = 1;
    
    std::string histogramFileName = argv[iFile];
    
    if (histogramFileName.find_first_of(':') != std::string::npos) {

      std::string aux = histogramFileName.substr(histogramFileName.find_first_of(':')+1);
      
      histogramFileName = histogramFileName.substr(0,histogramFileName.find_first_of(':'));
      
      if (aux.find_first_of(':') != std::string::npos) {
        maxShowers = std::stof(aux.substr(0,aux.find_first_of(':')));
        weight = std::stof( aux.substr(aux.find_first_of(':')+1) );
      } else {
        maxShowers = std::stof(aux);
      }

    }
    
    std::cout << "+ Parsing " << maxShowers << " showers in file " << histogramFileName << " (weight = " << weight << ")" << std::endl;

    CerenkovHistogramFile histogramFile(histogramFileName);
    
    if (histogramFile.IsZombie()) {
      continue; 
    }
      
    int numberOfShowersRead = 0;
    
    while(histogramFile.NextShower() && numberOfShowersRead++ < maxShowers) {
      while(histogramFile.NextHistogram()) {
        fitter.AddHistogram(histogramFile, 1.0/histogramFile.GetShowerEnergy());
      }
    }

  }
  
  fitter.Fit(parameterSeed);
  
  return 0;
}
