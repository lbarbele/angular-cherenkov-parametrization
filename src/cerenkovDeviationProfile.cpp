#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>

#include <TH2.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TCanvas.h>

#include <cerenkovDeviationProfile.h>


CerenkovDeviationProfile::CerenkovDeviationProfile(std::string deviationProfileName, int npt, double min, double max) :
  name(deviationProfileName),
  numberOfPoints(npt),
  minTheta(min),
  maxTheta(max),
  maxThetaPlot(max),
{
  // Nothing to do here.
}

CerenkovDeviationProfile::~CerenkovDeviationProfile()
{
  // Nothing to do here.
}

void CerenkovDeviationProfile::AddNewProfile(double age)
{
  int roundAge = std::round(age*10.0);
  
  if (this->mapOfProfiles.find(roundAge) == this->mapOfProfiles.end()) {

    this->mapOfProfiles.emplace(
      std::piecewise_construct ,
      std::forward_as_tuple(roundAge) ,
      std::forward_as_tuple((this->name + std::to_string(roundAge)).c_str(),
        "",
        this->numberOfPoints,
        this->minTheta,
        this->maxTheta,
        "s" )
    );

  }
}

int CerenkovDeviationProfile::CheckAge(const CerenkovHistogramFile & histogramFile)
{
  double ageLeft  = 10 * (histogramFile.GetAge() - histogramFile.GetDeltaAge());
  double ageRight = 10 * (histogramFile.GetAge() + histogramFile.GetDeltaAge());
  
  int roundAge;
  for (auto & pairAgeProfile : this->mapOfProfiles) {
    roundAge = pairAgeProfile.first;
    if (ageLeft < roundAge && roundAge < ageRight) {
      return roundAge;
    }
  }
  
  return -1;
}

void CerenkovDeviationProfile::AddData(const CerenkovHistogramFile & histogramFile, const CerenkovAngularDistribution & angularDist)
{
  // Check age
  int roundAge = this->CheckAge(histogramFile);
  if (roundAge <= 0) {
    return;
  }

  // Get a copy of the data histogram
  std::unique_ptr<TH1> dataHistogram((TH1*) histogramFile.GetThetaHistogramPointer()->Clone("copy") );
  dataHistogram->Scale(1.0/dataHistogram->Integral("width"));
  
  // Build a histogram with the fit model (already normalized)
  std::unique_ptr<TH1> fitHistogram( angularDist.GetFitHistogram(600,0,std::acos(-1)/3.0,dataHistogram->Integral()) );
  
  // Check binning
  if (dataHistogram->GetNbinsX() != this->numberOfPoints || fitHistogram->GetNbinsX() != this->numberOfPoints) {
    std::cerr << "CerenkovDeviationProfile::AddData: check histogram binning and try again" << std::endl;
    return;
  }
  
  // Get reference to correct profile according to shower age
  auto & profile = this->mapOfProfiles[roundAge];
  
  // Loop over points to fill the profile
  for (int ipt = 1; ipt <= this->numberOfPoints; ipt++) {
    double theta = profile.GetBinCenter(ipt);
    double dev = 1.0 - dataHistogram->GetBinContent(ipt)/fitHistogram->GetBinContent(ipt);
    profile.Fill(theta, dev);
  }
}

// This version uses a TH1 that could be any function as fit model
void CerenkovDeviationProfile::AddData(const CerenkovHistogramFile & histogramFile, TH1 * inputFitHistogram)
{
  // Check if fit histogram exists
  if (!inputFitHistogram) {
    return;
  }
  
  // Check age
  int roundAge = this->CheckAge(histogramFile);
  if (roundAge <= 0) {
    return;
  }
  
  // Get a copy of the data histogram
  std::unique_ptr<TH1> dataHistogram((TH1*) histogramFile.GetThetaHistogramPointer()->Clone("copy") );
  dataHistogram->Scale(1.0/dataHistogram->Integral("width"));
  
  // Make a copy of the fit histogram (supposed to be already normalized)
  std::unique_ptr<TH1D> fitHistogram( new TH1D(*(TH1D*)inputFitHistogram) );
  
  // Check axes
  if (dataHistogram->GetNbinsX() != this->numberOfPoints || fitHistogram->GetNbinsX() != this->numberOfPoints) {
    std::cerr << "CerenkovDeviationProfile::AddData: check histogram binning and try again" << std::endl;
    return;
  }
  
  // Get reference to correct profile according to shower age
  auto & profile = this->mapOfProfiles[roundAge];
  
  // Loop over points to fill the profile
  for (int ipt = 1; ipt <= this->numberOfPoints; ipt++) {
    double theta = profile.GetBinCenter(ipt);
    double dev = 1.0 - dataHistogram->GetBinContent(ipt)/fitHistogram->GetBinContent(ipt);
    profile.Fill(theta, dev);
  }
  
  // done!
}

void CerenkovDeviationProfile::Print(std::string outputPath)
{
  gSystem->mkdir(outputPath.c_str(),true);
  
  std::string outputFileName = outputPath + "/" + this->name + ".pdf";
  
  TCanvas canvas("","",600,300);
  
  canvas.Print((outputFileName+"[").c_str());
  
  TH2D hAxis("haxis","",100,this->minTheta,this->maxThetaPlot,100,-0.5,0.5);
  hAxis.SetStats(false);
  
  TLatex latexAgeLabel;
  
  for (auto & pairAgeProfile : this->mapOfProfiles) {

    auto & roundAge = pairAgeProfile.first;
    auto & profile = pairAgeProfile.second;
    
    std::string stringOfAge = "s = " + std::to_string(roundAge/10) + "." + std::to_string(roundAge%10);
    
    canvas.Clear();
    
    hAxis.Draw("axis");
    profile.Draw("same hist l");
    
    latexAgeLabel.DrawLatexNDC(0.05,0.95,stringOfAge.c_str());
    
    canvas.Print(outputFileName.c_str());
  }

  canvas.Print((outputFileName+"]").c_str());
}

void CerenkovDeviationProfile::WriteTable(std::string outputPath, int rebinFactor)
{
  // Create the output directory, if necessary
  gSystem->mkdir(outputPath.c_str(),true);
  
  // Open the output stream
  std::ofstream outStream(outputPath + "/" + this->name + ".dat");
  
  // Print table header
  outStream << std::setw(13) << "theta";
  for (auto & ppp : this->mapOfProfiles) {
    outStream << std::setw(13) << "s_" + std::to_string(ppp.first/10) + "." + std::to_string(ppp.first%10);
  }
  outStream << std::endl;
  
  // Print contents
  if (rebinFactor > 1 && this->numberOfPoints%rebinFactor == 0) {

  	auto rebinedProfiles = this->mapOfProfiles;

  	for (auto & ppp : rebinedProfiles) {
  		ppp.second.Rebin(rebinFactor);
    }
  		
  	auto & referenceProfile = (*rebinedProfiles.begin()).second;
  		
  	for (int iBin = 1; iBin <= referenceProfile.GetNbinsX(); iBin++) {
		  outStream << std::setw(13) << referenceProfile.GetBinCenter(iBin)*57.295779513;
		  for (auto & ppp : rebinedProfiles) {
        outStream << std::setw(13) << ppp.second.GetBinContent(iBin);
      }
		  outStream << std::endl;
		}

  } else {

  	auto & referenceProfile = (*this->mapOfProfiles.begin()).second;
  
		for (int iBin = 1; iBin <= referenceProfile.GetNbinsX(); iBin++) {
		  outStream << std::setw(13) << referenceProfile.GetBinCenter(iBin)*57.295779513;
		  for (auto & ppp : this->mapOfProfiles) {
        outStream << std::setw(13) << ppp.second.GetBinContent(iBin);
      }
		  outStream << std::endl;
		}

  }
  
}
