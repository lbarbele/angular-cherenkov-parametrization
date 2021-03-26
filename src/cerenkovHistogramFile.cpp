#include <iostream>
#include <string>

#include <TTree.h>

#include <cerenkovHistogramFile.h>

/**
 * Constructor with given path
 */
CerenkovHistogramFile::CerenkovHistogramFile(std::string histogramFileName) :
  TFile(histogramFileName.c_str(),"read"),
  numberOfShowers(0),
  nextShower(0),
  eventTreePointer(nullptr),
  runTreePointer(nullptr),
  eventDirectoryPointer(nullptr),
  rhoHistogram(nullptr),
  ageHistogram(nullptr),
  refHistogram(nullptr),
  nextHistogram(0),
  numberOfHistograms(0),
  currentThetaHistogram(nullptr),
  currentRunNumber(0),
  currentEventNumber(0),
  currentAge(0),
  currentDeltaAge(0),
  currentRefractiveIndex(0),
  currentDensity(0),
  currentAtmosphericBin(0),
  showerEnergy(0),
  showerZenith(0),
  showerAzimuth(0)
{
  // Check if file was correctly opened
  if (this->IsZombie()) {
    std::cerr << "CerenkovHistogramFile::CerenkovHistogramFile: unable to open file " << histogramFileName << "." << std::endl;
    return;
  }
  
  // Try to get pointer the event tree
  this->eventTreePointer.reset( this->Get<TTree>("tEvt") );
  
  // Check pointer to the event tree
  if (!this->eventTreePointer) {
    std::cerr << "CerenkovHistogramFile::CerenkovHistogramFile: unable to read event tree from file" << histogramFileName << "." << std::endl;
    return;
  }
  
  // Try to get pointer to the run tree
  this->runTreePointer.reset( this->Get<TTree>("tRun") );
  
  // Check pointer to the run tree
  if (!this->runTreePointer) {
    std::cerr << "CerenkovHistogramFile::CerenkovHistogramFile: unable to read run tree from file" << histogramFileName << "." << std::endl;
    return;
  }

  // Read both branches (header/end) from the event tree
  this->eventTreePointer->SetBranchAddress("fEvtHeader",this->eventHeaderArray.data());
  this->eventTreePointer->SetBranchAddress("fEvtEnd",this->eventEndArray.data());
  
  // Read both branches (header/end) from the run tree 
  this->runTreePointer->SetBranchAddress("fRunHeader",this->runHeaderArray.data());
  this->runTreePointer->SetBranchAddress("fRunEnd",this->runEndArray.data());
  
  // Get number of showers in current file
  this->numberOfShowers = this->eventTreePointer->GetEntries();
}

/**
 * Read the next shower
 */
bool CerenkovHistogramFile::NextShower()
{
  // Reset all data related to last shower read
  this->nextHistogram = 0;
  this->numberOfHistograms = 0;
  this->showerEnergy = 0;
  this->showerZenith = 0;
  this->currentAtmosphericBin = 0;
  this->rhoHistogram.reset(nullptr);
  this->ageHistogram.reset(nullptr);
  this->refHistogram.reset(nullptr);
  this->currentThetaHistogram.reset(nullptr);
  this->eventDirectoryPointer.reset(nullptr);
  
  // Check if there is a next shower. If there is not, return false
  if (this->nextShower >= this->numberOfShowers) {
    return false;
  }
    
  // Read the next shower in event tree
  this->eventTreePointer->GetEntry(this->nextShower);
  
  // Get run/event number
  this->currentRunNumber = this->eventHeaderArray[43];
  this->currentEventNumber = this->eventHeaderArray[1];
  
  // Get pointer to event directory
  std::string eventDirectoryName = "Event_" + std::to_string(this->currentRunNumber) + "_" + std::to_string(this->currentEventNumber);
  this->eventDirectoryPointer.reset( this->GetDirectory(eventDirectoryName.c_str()) );
  
  // Check if the pointer is not null
  if (!this->eventDirectoryPointer) {
    return false;
  }

  // Get shower information
  this->showerEnergy  = this->eventHeaderArray[3]*1e-3;
  this->showerZenith  = this->eventHeaderArray[10];
  this->showerAzimuth = this->eventHeaderArray[11];
    
  // Get pointers to density, refractive index and age (hisogram) profiles
  this->rhoHistogram.reset( (TH1*) this->eventDirectoryPointer->Get("Profiles/hDensity") );
  this->ageHistogram.reset( (TH1*) this->eventDirectoryPointer->Get("Profiles/hAge") );
  this->refHistogram.reset( (TH1*) this->eventDirectoryPointer->Get("Profiles/hReffrac") );
  
  // Set next histogram index and number of histograms
  this->numberOfHistograms = this->ageHistogram->GetNbinsX();
  this->nextHistogram = 1;
  
  // Increment nextShower counter
  this->nextShower++;
  
  return true;
}

/**
 * Read shower with a given index
 */
bool CerenkovHistogramFile::ReadShower(int showerId)
{
  this->nextShower = showerId;
  return this->NextShower();
}

/**
 * Read next histogram of current shower
 */
bool CerenkovHistogramFile::NextHistogram()
{
  // Reset previous histogram
  this->currentThetaHistogram.reset(nullptr);
  this->currentRefractiveIndex = 0;
  this->currentAge = 0;
  this->currentDeltaAge = 0;
  this->currentDensity = 0;
  this->currentAtmosphericBin = 0;
  
  // Check if there is a next histogram
  if (this->nextHistogram > this->numberOfHistograms) {
    return false;
  }
    
  // Get next histogram
  this->currentThetaHistogram.reset( (TH1*) this->eventDirectoryPointer->Get(("hTheta_"+std::to_string(nextHistogram)).c_str()) );
  
  // Check if pointer to histogram is not null
  if (!this->currentThetaHistogram) {
    return false;
  }
    
  // Get age/refractive index
  this->currentRefractiveIndex = this->refHistogram->GetBinContent(this->nextHistogram);
  this->currentAge = this->ageHistogram->GetBinContent(this->nextHistogram);
  this->currentDeltaAge = this->ageHistogram->GetBinError(this->nextHistogram);
  this->currentDensity = this->rhoHistogram->GetBinContent(this->nextHistogram);
  this->currentAtmosphericBin = this->nextHistogram;
  
  // Increment histogram counter
  this->nextHistogram++;

  return true;
}

/**
 * Read histogram of given atmospheric bin
 */
bool CerenkovHistogramFile::ReadHistogram(int histogramId)
{
  this->nextHistogram = histogramId;
  return this->NextHistogram();
}




/**
 * Initiliaze the static arrays with information of the atmospheric model
 */
const double CerenkovHistogramFile::vAtmH[5] = {  -577950,  400000,   1e+06,   4e+06,     1e+07};

const double CerenkovHistogramFile::vAtmA[5] = { -186.556, -94.919, 0.61289,       0, 0.0112829};

const double CerenkovHistogramFile::vAtmB[5] = {  1222.66, 1144.91, 1305.59, 540.178,         1};

const double CerenkovHistogramFile::vAtmC[5] = {   994186,  878154,  636143,  772170,     1e+09};

const double CerenkovHistogramFile::vMeanH[20] = { 27.4575, 18.3503, 14.9991,  12.8283, 11.2159
                                                 , 9.92651, 8.81325, 7.82429,  6.93561, 6.12871
                                                 , 5.38979, 4.70828, 4.07563,  3.48148, 2.92049
                                                 , 2.38947, 1.88538, 1.40563, 0.947979, 0.51047 };

/**
 * Compute atmospheric depth
 */
double CerenkovHistogramFile::AtmDepth(double height)
{
  for (int iLayer = 0; iLayer < 4; iLayer++) {
    if (height < vAtmH[iLayer+1]) {
      return vAtmA[iLayer] + vAtmB[iLayer] * std::exp( - height / vAtmC[iLayer] );
    }
  }
  return std::max(0.0, vAtmA[4] - (vAtmB[4]/vAtmC[4]) * height);
}

/**
 * Compute atmospheric height
 */
double CerenkovHistogramFile::AtmHeight(double depth)
{
  for (int iLayer = 0; iLayer < 4; iLayer++) {
    if (depth > AtmDepth(vAtmH[iLayer+1])) {
      return vAtmC[iLayer] * std::log(vAtmB[iLayer]/(depth - vAtmA[iLayer]));
    }
  }
  return (vAtmC[4]/vAtmB[4]) * (vAtmA[4] - depth);
}

/**
 * Compute refractive index minus 1
 */
double CerenkovHistogramFile::AtmDelta(double height)
{
  for (int iLayer = 0; iLayer < 4; iLayer++) {
    if (height < vAtmH[iLayer+1]) {
      return 0.000283*std::exp( - height / vAtmC[iLayer] );
    }
  }
  return 0.000283*std::exp( - height / vAtmC[4] );
}

/**
 * Compute atmospheric density
 */
double CerenkovHistogramFile::AtmRho(double height)
{
  for (int iLayer = 0; iLayer < 4; iLayer++) {
    if (height < vAtmH[iLayer+1]) {
      return (vAtmB[iLayer]/vAtmC[iLayer]) * std::exp( - height / vAtmC[iLayer] );
    }
  }
  return std::max(0.0,  (vAtmB[4]/vAtmC[4]) * height);
}
