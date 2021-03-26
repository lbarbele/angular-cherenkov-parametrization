#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include <TSystem.h>

#include <cerenkovCanvas.h>
#include <cerenkovHistogramFile.h>
#include <cerenkovDeviationProfile.h>
#include <cerenkovAngularDistribution.h>

int main(int argc, char ** argv)
{
  int numberOfThetaBins = 600;
  double minTheta = 0;
  double maxTheta = std::acos(-1)/3.0;
  int deviationRebinFactor = 6;
  
  int maxPlotShowers = 10;
  
  std::vector<double> ageVector = {0.8,1.0,1.2};
  
  // Set the parameter vector - these are the final results!!
  std::vector<double> p_gamma  = {0.34329, -0.10683, 1.46852, 1.4053, 0.32382, -0.048841, 1.22206, 0.26472, 0.0031206, 0.       , 0.     , 0.      };
  std::vector<double> p_proton = {0.21155, -0.16639, 1.21803, 4.513 , 0.45092, -0.058687, 1.32447, 0.41722, 0.009528 , 0.022552, -0.4207, -0.008843};
  
  // get the output prefix
  std::string outputPrefix = argv[1];
  std::string plotDirectory = "plots/test-parametrization/" + (deviationRebinFactor == 1 ? outputPrefix : "rebin_" + outputPrefix);
  std::string tableDirectory = "tables/test-parametrization/" + (deviationRebinFactor == 1 ? outputPrefix : "rebin_" + outputPrefix);
  gSystem->mkdir(plotDirectory.c_str(),true);
  gSystem->mkdir(tableDirectory.c_str(),true);
  
  // Create the deviation profile using the first input parameter as name
  CerenkovDeviationProfile deviationProfile("deviationProfile",numberOfThetaBins,minTheta,maxTheta);
  CerenkovDeviationProfile lemoineProfile("lemoineProfile",numberOfThetaBins,minTheta,maxTheta);
  CerenkovDeviationProfile gillerProfile("gillerProfile",numberOfThetaBins,minTheta,maxTheta);
  CerenkovDeviationProfile nerlingProfile("nerlingProfile",numberOfThetaBins,minTheta,maxTheta);
  
  // Set ages to be plotted
  for (auto age : ageVector) {
    deviationProfile.AddNewProfile(age);
    lemoineProfile.AddNewProfile(age);
    gillerProfile.AddNewProfile(age);
    nerlingProfile.AddNewProfile(age);
  }
  
  // Set shower counter
  int iShower = 0;
  
  // loop over input files
  for (int iFile = 2; iFile < argc; iFile++) {

    // open the histogram file
    std::string histogramFileName = argv[iFile];
    CerenkovHistogramFile histogramFile(histogramFileName);
    
    // check file
    if (histogramFile.IsZombie()) {
      std::cerr << "Skipping file " << histogramFileName << std::endl;
      continue;
    }
    
    // Set parameter vector according to primary
    std::vector<double> p;
    
    if (histogramFileName.find("_gamma_") != std::string::npos) {
      p = p_gamma;
    } else if (histogramFileName.find("_p_") != std::string::npos) {
      p = p_proton;
    } else {
      std::cerr << "Unable to guess primary particle from file " << histogramFileName << ". It will be skipped!" << std::endl;
      continue;
    }
    
    // loop over showers
    while(histogramFile.NextShower()) {
     
      // cut bad showers
      if ( ( histogramFile.GetRunNumber() == 9 )
           && ( histogramFile.GetEvtNumber() == 4 )
           && ( histogramFileName == "data/histograms/vertical_veritas_QGSII_urqmd_gamma_1TeV.root" ) 
         ) {
        continue;
      }

      // Create the canvas for this shower only if maxPlotShowers has not been reached
      std::string plotFileName = plotDirectory + "/shower_" + std::to_string(iShower+1) + ".pdf";
      std::unique_ptr<CerenkovCanvas> showerCanvas( iShower < maxPlotShowers ? new CerenkovCanvas(plotFileName) : nullptr );

      // loop over histograms of this shower
      while(histogramFile.NextHistogram()) {
       
        // read data of this atmospheric bin
        double refractiveIndex = histogramFile.GetRefractiveIndex();
        double delta = refractiveIndex - 1;
        double emissionAngle = std::acos(1.0/refractiveIndex);
        double age = histogramFile.GetAge();
        double deltaAge = histogramFile.GetDeltaAge();
        double energy = histogramFile.GetShowerEnergy();

        // make age string and skip ages we don't want to plot
        std::string ageString = "";

        for (auto referenceAge : ageVector) {
          if (std::fabs(age-referenceAge) <= deltaAge) {
            ageString = "s_" + std::to_string(referenceAge).substr(0,3);
          }
        }

        if (ageString == "") {
          continue;
        }

        // Compute fit parameters
        double parNu     = p[0] * std::pow(delta, p[1]) + p[2] * std::log(age);
        double parTheta1 = p[3] * std::pow(delta, p[4]) * std::pow(energy,p[11]) + p[5] * std::log(age);
        double parRatio  = p[6] + p[7]*(age-1);
        double parEps    = p[8] + p[9]*std::pow(energy,p[10]);
        
        // Create the angular distribution from computed parameters
        CerenkovAngularDistribution angularDist;
        angularDist.SetParameters(emissionAngle, parTheta1, parRatio, parNu, parEps);
        
        // Add this angular distribution to the deviation profile
        deviationProfile.AddData(histogramFile, angularDist);
        
        // Add data to deviation profiles from reference parametrizations
        std::unique_ptr<TH1> lemoineHistogram( CerenkovAngularDistribution::GetLemoineHistogram(histogramFile, numberOfThetaBins) );
        std::unique_ptr<TH1> gillerHistogram( CerenkovAngularDistribution::GetGillerHistogram(histogramFile, numberOfThetaBins) );
        std::unique_ptr<TH1> nerlingHistogram( CerenkovAngularDistribution::GetNerlingHistogram(histogramFile, numberOfThetaBins) );
        lemoineProfile.AddData(histogramFile, lemoineHistogram.get());
        gillerProfile.AddData(histogramFile, gillerHistogram.get());
        nerlingProfile.AddData(histogramFile, nerlingHistogram.get());
        
        // Add this to shower canvas, if it exists
        if (showerCanvas) {
          showerCanvas->AddPage(&histogramFile, &angularDist);
          showerCanvas->MakeDataTable("data" ,ageString,&histogramFile);
          showerCanvas->MakeDataTable("fit-a",ageString,&angularDist,4000);
          showerCanvas->MakeDataTable("fit-b",ageString,&angularDist, 250);
          showerCanvas->MakeDataTable("lemoine",ageString,lemoineHistogram.get());
          showerCanvas->MakeDataTable("giller",ageString,gillerHistogram.get());
          showerCanvas->MakeDataTable("nerling",ageString,nerlingHistogram.get());
        }
      }
      
      // Print tables, if existent
      if (showerCanvas) {
        showerCanvas->PrintTable("data" ,std::ofstream( tableDirectory + "/shower_" + std::to_string(iShower+1) + "_data.dat" ));
        showerCanvas->PrintTable("fit-a",std::ofstream( tableDirectory + "/shower_" + std::to_string(iShower+1) + "_fit-a.dat" ), 8);
        showerCanvas->PrintTable("fit-b",std::ofstream( tableDirectory + "/shower_" + std::to_string(iShower+1) + "_fit-b.dat" ));
        showerCanvas->PrintTable("lemoine",std::ofstream( tableDirectory + "/shower_" + std::to_string(iShower+1) + "_lemoine.dat" ));
        showerCanvas->PrintTable("giller",std::ofstream( tableDirectory + "/shower_" + std::to_string(iShower+1) + "_giller.dat" ));
        showerCanvas->PrintTable("nerling",std::ofstream( tableDirectory + "/shower_" + std::to_string(iShower+1) + "_nerling.dat" ));
      }
      
      // Increment shower counter
      iShower++;
    }
  }
    
  // Print the deviation profiles
  deviationProfile.Print(plotDirectory);
  lemoineProfile.Print(plotDirectory);
  gillerProfile.Print(plotDirectory);
  nerlingProfile.Print(plotDirectory);
  
  // Write tables of deviation profiles
  deviationProfile.WriteTable(tableDirectory,deviationRebinFactor);
  lemoineProfile.WriteTable(tableDirectory,deviationRebinFactor);
  gillerProfile.WriteTable(tableDirectory,deviationRebinFactor);
  nerlingProfile.WriteTable(tableDirectory,deviationRebinFactor);

  return 0;
}
