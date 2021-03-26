#include <iostream>
#include <iomanip>
#include <future>
#include <cmath>

#include <TH1.h>
#include <Fit/Fitter.h>
#include <Fit/FitResult.h>

#include <gsl/gsl_sf_gamma.h>

#include <cerenkovFitter.h>

CerenkovFitter::CerenkovFitter() :
  headOfEntryList(nullptr),
  tailOfEntryList(nullptr),
  nextEntryFromList(nullptr),
  parCalcPtr(nullptr),
  minAge(0),
  maxAge(3),
  isVerbose(false)
{
  // Nothing to do here
}

// Note that the pointer to the next histogram has to be set externally.
bool CerenkovFitter::AddHistogram(const CerenkovHistogramFile & cerHistogramFileRef, double normalizationFactor)
{
  // Check if age is inside specified interval
  double age = cerHistogramFileRef.GetAge();
  double deltaAge = cerHistogramFileRef.GetDeltaAge();
  
  if (age+deltaAge < this->minAge || age-deltaAge > this->maxAge) {
    return false;
  }

  // Make a copy of the histogram
  TH1D h(*(TH1D*)cerHistogramFileRef.GetThetaHistogramPointer());
  
  // Normalize the histogram, if asked for
  if (normalizationFactor != 1 && normalizationFactor > 0) {
    h.Scale(normalizationFactor);
  }
  
  // Create a new data entry on the list and assign it as the tail of the list
  this->tailOfEntryList = new CerenkovFitter::DataEntry( cerHistogramFileRef.GetShowerEnergy()
                                                       , cerHistogramFileRef.GetShowerZenith()
                                                       , cerHistogramFileRef.GetAge()
                                                       , cerHistogramFileRef.GetRefractiveIndex()
                                                       , h
                                                       , this->tailOfEntryList );
  
  // If the list has no head, this is the first entry and the head corresponds to the tail
  if (!this->headOfEntryList) {
    this->headOfEntryList = this->tailOfEntryList;
  }
  
  return true;
}

ROOT::Fit::FitResult CerenkovFitter::Fit(std::vector<double> seedVector, std::vector<double> lowerVector, std::vector<double> upperVector)
{
  // Create the root fitter
  ROOT::Fit::Fitter rootFitter;
  
  this->numberOfParameters = seedVector.size();
  
  // Create parameter configuration and copy given seeds
  rootFitter.Config().SetParamsSettings(seedVector.size(),seedVector.data());
  
  // Check if bounds are given and set them - if lower >= up, the parameter will be fixed
  if (lowerVector.size() == seedVector.size() && upperVector.size() == seedVector.size()) {
    for (int ipar = 0; ipar < seedVector.size(); ipar++) {
      if (lowerVector[ipar] >= upperVector[ipar]) {
        rootFitter.Config().ParSettings(ipar).Fix();
      } else {
        rootFitter.Config().ParSettings(ipar).SetLimits(lowerVector[ipar],upperVector[ipar]);
      }
    }
  }
  
  // Set objective function
  rootFitter.SetFCN<CerenkovFitter>(seedVector.size(),*this);

  // Perform the fit
  rootFitter.FitFCN();
  
  // Print parameters, if verbosity flag is raised
  if (this->isVerbose) {
    rootFitter.Result().Print(std::cout);
    for (int ipar = 0; ipar < seedVector.size(); ipar++) {
      std::cout << rootFitter.Result().Parameter(ipar) << ", ";
    }
    std::cout << std::endl;
  }
  
  return rootFitter.Result();
}

double CerenkovFitter::operator()(const double * p)
{
  // Set the first entry from the list as the next one
  this->nextEntryFromList = this->headOfEntryList;
  
  // This is the vector of future<doubles> that will hold partial results computed at each thread
  std::vector<std::future<double>> vectorOfFutureDeviances;
  
  // Launc nproc threads to compute deviance in parallel
  int nproc = 4;
  for (int i=0; i<nproc; i++) {
    vectorOfFutureDeviances.emplace_back( std::async(&CerenkovFitter::ComputeDevianceAsync,this,p) );
  }
    
  // Retrieve results and sum all contributions
  double deviance = 0;
  for (auto & futureDeviance : vectorOfFutureDeviances) {
    deviance += futureDeviance.get();
  }
  
  // Print, if verbosity is enabled
  if (this->isVerbose) {
    std::cout << std::setw(13) << deviance;
    for (int ipar = 0; ipar < this->numberOfParameters; ipar++) {
      std::cout << std::setw(13) << p[ipar];
    }
    std::cout << std::endl;
  }

  // Return deviance value
  return deviance;
}

CerenkovFitter::DataEntry * CerenkovFitter::GetNextEntryFromListAsync()
{
  // Lock the mutex and block other threads
  this->nextEntryReaderBarrier.lock();
  
  // This is the actual next entry (may be a null pointer)
  auto returnValue = this->nextEntryFromList;
  
  // Set the new next entry, if current next exists
  if (this->nextEntryFromList) {
    this->nextEntryFromList = this->nextEntryFromList->nextEntry;
  }
    
  // Unlock the mutex, allowing other threads to read the next entry
  this->nextEntryReaderBarrier.unlock();
  
  // Return the next entry
  return returnValue;
}



double CerenkovFitter::ComputeDevianceAsync(const double * parameterArrayPointer)
{
  double deviance = 0;
  
  double l1, l2, lastTerm, nuMinusOne;
  
  double functionParameterArray[4];
  double & t1  = functionParameterArray[0];
  double & t2  = functionParameterArray[1];
  double & nu  = functionParameterArray[2];
  double & eps = functionParameterArray[3];
  
  while(auto currentEntry = this->GetNextEntryFromListAsync()) {
    double & age = currentEntry->showerAge;
    
    this->parCalcPtr(functionParameterArray, parameterArrayPointer, currentEntry->showerAge, currentEntry->refractiveIndex-1, currentEntry->showerEnergy);
    
    l1 = 1.0/t1;
    l2 = 1.0/t2;
    nuMinusOne = nu-1;
    
    double M = 0;
    
    for (int k=0; k < currentEntry->numberOfBins; k++) {
      lastTerm = 1 + eps * std::exp(currentEntry->vector_x[k]*l2);
      deviance -= currentEntry->vector_y[k] * std::log( lastTerm );
      M += currentEntry->vector_j[k] * std::pow(currentEntry->vector_x[k],nuMinusOne) * std::exp(-currentEntry->vector_x[k]*l1) * lastTerm;
    }
    
    deviance += currentEntry->N * std::log(M / currentEntry->N);
    deviance += currentEntry->A;
    deviance -= currentEntry->B * nuMinusOne;
    deviance += currentEntry->C * l1;
  }
  
  deviance *= 2;
  return deviance;
}

//
// The CerenkovFitter::DatEntry class
//

// Constructor
CerenkovFitter::DataEntry::DataEntry( double showerEnergyInput
                                    , double showerZenithInput
                                    , double showerAgeInput
                                    , double refractiveIndexInput
                                    , const TH1 & histogram, CerenkovFitter::DataEntry * previousEntryPointerInput) :
  showerEnergy(showerEnergyInput),
  showerZenith(showerZenithInput),
  showerAge(showerAgeInput),
  refractiveIndex(refractiveIndexInput),
  emissionAngle(std::acos(1.0/refractiveIndexInput)),
  minTheta(histogram.GetBinLowEdge(1)),
  maxTheta(histogram.GetBinLowEdge(histogram.GetNbinsX()+1)),
  numberOfBins(histogram.GetNbinsX()),
  vector_y(histogram.GetNbinsX(),0),
  vector_x(histogram.GetNbinsX(),0),
  vector_j(histogram.GetNbinsX(),0),
  A(0),
  B(0),
  C(0),
  N(0),
  previousEntry(previousEntryPointerInput),
  nextEntry(nullptr)
{
  static const double pi = std::acos(-1);
  
  // Loop over bins of input histogram to get v and y vectors
  for (int iBin = 1; iBin <= histogram.GetNbinsX(); iBin++) {

    // compute <theta_p> = max(theta,tch)
    double theta = histogram.GetBinCenter(iBin);
    double deltaTheta = 0.5 * histogram.GetBinWidth(iBin);
    double thetap = std::max(this->emissionAngle, theta);
    double entries = histogram.GetBinContent(iBin);
    
    // The peak function is the average of (sin(theta)/sin(theta_p)) * ( pi - log(1 - min(theta,theta_p)/max(theta,theta_p)) )
    double xk = theta/this->emissionAngle;
    double lk = deltaTheta/this->emissionAngle;
    double peakFunction = 0;
    if (xk+lk <= 1) {
      peakFunction += xk*lk * (2*pi + 1 - std::log(1+xk*xk-2*xk-lk*lk) ) + lk + (xk*xk+lk*lk-1) * std::atanh(lk/(1-xk));
    } else if (xk-lk < 1) {
      peakFunction += pi + 1.5 - (xk*xk + lk*lk - 2*xk*lk)*(pi+0.5) - (xk-lk) - (1-xk*xk-lk*lk+2*xk*lk)*std::log(1-xk+lk);
      peakFunction *= 0.5;
      peakFunction += (xk+lk-1)*(pi - std::log(xk+lk-1)) + (xk+lk)*std::log(xk+lk);
    } else {
      peakFunction += lk*(pi+std::log(std::sqrt( (xk*xk-lk*lk) / (xk*xk-2*xk+1-lk*lk) ))) + xk*std::atanh(lk/xk) - (xk-1)*std::atanh(lk/(xk-1));
      peakFunction *= 2;
    } 
    peakFunction *= this->emissionAngle;
    
    // x is theta_p
    this->vector_x[iBin-1] = thetap;
    
    // j is peakFunction
    this->vector_j[iBin-1] = peakFunction;
    
    // y is the vector of bin contents
    this->vector_y[iBin-1] = entries;
    
    if (entries == 0) continue;
    
    // N is the sum of entries
    this->N += entries;
    
    // A, B, and C are more complex functions of data
    this->A += entries * std::log(entries/peakFunction);
    this->B += entries * std::log(thetap);
    this->C += entries * thetap;
  }
  
  // If there is a previous entry, tell it this entry is the next
  if (previousEntryPointerInput) {
    previousEntryPointerInput->nextEntry = this;
  }
}



// The destructor must ensure that if an element in the middle of the list is destructed,
// the remaining elements are still connected
CerenkovFitter::DataEntry::~DataEntry()
{
  // If there is a previous entry, the new nextEntry of the previousEntry is the nextEntry of this object
  if (this->previousEntry) {
    this->previousEntry->nextEntry = this->nextEntry;
  }
  
  // Likewise, if there is a nextEntry, the new previousEntry of the nextEntry is the previousEntry of this object
  if (this->nextEntry) {
    this->nextEntry->previousEntry = this->previousEntry;
  }
    
  // The other fields can be default destructed
}
