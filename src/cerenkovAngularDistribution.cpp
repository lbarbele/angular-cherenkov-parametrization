#include <cmath>
#include <vector>
#include <iomanip>
#include <iostream>

#include <TH1.h>
#include <Fit/FitResult.h>

#include <cerenkovAngularDistribution.h>
#include <cerenkovFitter.h>


// Constructor: do nothing
CerenkovAngularDistribution::CerenkovAngularDistribution() :
  emissionAngle(0.02),
  parTheta_1(0.08),
  parRatio(1.2),
  parNu(0.8),
  parEps(0.0005),
  errTheta_1(0),
  errRatio(0),
  errNu(0),
  errEps(0),
  parameterSeedArray({parTheta_1,parRatio,parNu,parEps}),
  lowerParameterBound({0.00001, -10,  0, 0.00001}),
  upperParameterBound({      1,  10, 10,       1}),
{

}

void CerenkovAngularDistribution::SetParameterSeed(double t1, double r, double nu, double eps)
{
  this->parameterSeedArray = {t1,r,nu,eps};
}

void CerenkovAngularDistribution::SetLowerBounds(double t1, double r, double nu, double eps)
{
  this->lowerParameterBound = {t1,r,nu,eps};
}

void CerenkovAngularDistribution::SetUpperBounds(double t1, double r, double nu, double eps)
{
  this->upperParameterBound = {t1,r,nu,eps};
}

// Fit function: return a FitResult object
ROOT::Fit::FitResult CerenkovAngularDistribution::EstimateParameters(const CerenkovHistogramFile & cerHistogramFileRef)
{
  // Get emission angle
  this->emissionAngle = std::acos(1.0/cerHistogramFileRef.GetRefractiveIndex());
  
  // Build fitter and add data from histogram
  CerenkovFitter fitter;
  fitter.SetParameterFunction(CerenkovAngularDistribution::ParameterCalculator);
  fitter.AddHistogram(cerHistogramFileRef, false);
  
  // Do the fit and retrieve the result
  this->fitResult = fitter.Fit(this->parameterSeedArray, this->lowerParameterBound, this->upperParameterBound);
  
  // update parameters
  this->parTheta_1 = this->fitResult.Parameter(0);
  this->parRatio   = this->fitResult.Parameter(1);
  this->parNu      = this->fitResult.Parameter(2);
  this->parEps     = this->fitResult.Parameter(3);
  
  // update parameter errors
  this->errTheta_1 = this->fitResult.Error(0);
  this->errRatio   = this->fitResult.Error(1);
  this->errNu      = this->fitResult.Error(2);
  this->errEps     = this->fitResult.Error(3);
  
  // Return the fitresult object
  return this->fitResult;
}

// Set parameter (used as guess in the fit)
void CerenkovAngularDistribution::SetParameters(double tch, double t1, double r, double nu, double eps)
{
  this->emissionAngle = tch;
  this->parTheta_1    = t1;
  this->parRatio      = r;
  this->parNu         = nu;
  this->parEps        = eps;
}

// Evaluate a vector of numberOfPoints points in the range minTheta,maxTheta (return y values only)
std::vector<double> CerenkovAngularDistribution::EvalVector(int numberOfPoints, double minTheta, double maxTheta, double normalization ) const 
{
  static const double pi = std::acos(-1);
  
  std::vector<double> retVector(numberOfPoints,0);
  
  double dTheta = (maxTheta-minTheta)/double(numberOfPoints);
  double integral = 0;
  
  for (int iPoint = 0; iPoint < numberOfPoints; iPoint++) {
    
    double theta = minTheta + (0.5 + iPoint) * dTheta;
    double theta_p = std::max(theta, this->emissionAngle);
    
    double xk = theta/this->emissionAngle;
    double lk = 0.5*dTheta/this->emissionAngle;

    double peakFunction = 0;

    if (xk+lk <= 1) {
      peakFunction += xk*lk * (2*pi + 1 - std::log(1+xk*xk-2*xk-lk*lk) ) + lk + (xk*xk+lk*lk-1) * std::atanh(lk/(1-xk));
    } else if (xk-lk < 1)  {
      peakFunction += pi + 1.5 - (xk*xk + lk*lk - 2*xk*lk)*(pi+0.5) - (xk-lk) - (1-xk*xk-lk*lk+2*xk*lk)*std::log(1-xk+lk);
      peakFunction *= 0.5;
      peakFunction += (xk+lk-1)*(pi - std::log(xk+lk-1)) + (xk+lk)*std::log(xk+lk);
    } else {
      peakFunction += lk*(pi+std::log(std::sqrt( (xk*xk-lk*lk) / (xk*xk-2*xk+1-lk*lk) ))) + xk*std::atanh(lk/xk) - (xk-1)*std::atanh(lk/(xk-1));
      peakFunction *= 2;
    } 

    peakFunction *= this->emissionAngle;
    
    double functionValue = peakFunction;
    functionValue *= std::pow(theta_p,this->parNu-1);
    functionValue *= std::exp(-theta_p/this->parTheta_1);
    functionValue *= 1 + this->parEps * std::exp(theta_p/(this->parRatio*this->parTheta_1));
    
    integral += functionValue;
    
    retVector[iPoint] = functionValue;
  }
  
  for (int iPoint = 0; iPoint < numberOfPoints; iPoint++) {
    retVector[iPoint] *= normalization / integral;
  }
  
  return retVector;
}



// Create a histogram from fit function
TH1 * CerenkovAngularDistribution::GetFitHistogram(int numberOfPoints, double minTheta, double maxTheta, double normalization ) const
{
  auto yVector = this->EvalVector(numberOfPoints,minTheta,maxTheta,normalization);
  TH1 * hFit = new TH1D("","",numberOfPoints,minTheta,maxTheta);

  for (int ipt = 0; ipt < numberOfPoints; ipt++) {
    hFit->SetBinContent(ipt+1, yVector[ipt]);
  }

  return hFit;
}



//
// Static method
//
void CerenkovAngularDistribution::ParameterCalculator(double * functionParameters, const double *& p, const double &age, const double &delta, const double &energy)
{
  functionParameters[0] = p[0]; // t1 
  functionParameters[1] = functionParameters[0] * p[1]; // t2
  functionParameters[2] = p[2]; // nu
  functionParameters[3] = p[3]; // eps
}



//
//
//
TH1 * CerenkovAngularDistribution::GetLemoineHistogram(double showerZenith, int numberOfPoints)
{
  TH1 * hist = new TH1D("","",numberOfPoints,0,std::acos(-1)/3.0);
  
  double eta = 0.015 * std::sqrt(std::cos(showerZenith));
  
  for (int iBin = 1; iBin <= numberOfPoints; iBin++) {
   
    double theta = hist->GetBinCenter(iBin);
    double dTheta = 0.5*hist->GetBinWidth(iBin);
    
    // value at bin center
    hist->SetBinContent(iBin, (2.0/9.0) * ( theta < eta ? theta/eta : std::exp(- 0.25 * (theta/eta-1) ) ) / eta);
  }
  
  return hist;
}

TH1 * CerenkovAngularDistribution::GetGillerHistogram(double height, double age, int numberOfPoints)
{
  static const double pi = std::acos(-1);
  
  TH1 * hist = new TH1D("","",numberOfPoints,0,std::acos(-1)/3.0);
  
  double thetaZero =   6.058 * std::exp(-0.001103*age     - 0.002886*height) - 5.447;
  double a1        =   2.905 * std::exp( -0.03851*age     +   0.1072*height) + 10.66;
  double c1        =   7.320 * std::exp(  -0.3778*age     +  0.07202*height) + 7.143;
  double c2        = - 2.754 * std::exp(  -0.4242*age*age +   0.1084*height) + 3.344*age*age - 4.294*age;
  double alpha     =   2.620 * std::exp(  0.06837*age     -  0.02247*height) + 1.110;
  double a2        = a1 * std::pow(thetaZero,alpha) * std::exp( - c1*thetaZero - c2*thetaZero*thetaZero );
  
  for (int iBin = 1; iBin <= numberOfPoints; iBin++) {

    double theta = hist->GetBinCenter(iBin);
    double dTheta = 0.5*hist->GetBinWidth(iBin);

    // Use value at bin center
    hist->SetBinContent(iBin, (theta < thetaZero) ? a1 * std::exp(-c1*theta - c2*theta*theta) : a2 * std::pow(theta, -alpha));
  }
  
  return hist;
}

TH1 * CerenkovAngularDistribution::GetNerlingHistogram(double age, double thresholdEnergy, int numberOfPoints)
{
  TH1 * hist = new TH1D("","",numberOfPoints,0,std::acos(-1)/3.0);
  
  double a   =  0.42489 +  0.58371*age - 0.082373*age*age;
  double b   = 0.055108 - 0.095587*age + 0.056952*age*age;
  double tc  = 0.62694 * std::pow(thresholdEnergy, -0.60590);
  double tcc = tc * (10.509 - 4.9644*age);
  
  for (int iBin = 1; iBin <= numberOfPoints; iBin++) {
   
    double theta = hist->GetBinCenter(iBin);
    double dTheta = 0.5*hist->GetBinWidth(iBin);
    
    // value at bin center
    hist->SetBinContent(iBin, (a/tc)*std::exp(-theta/tc) + (b/tcc)*std::exp(-theta/tcc) );
  }
  
  return hist;
}

TH1 * CerenkovAngularDistribution::GetLemoineHistogram(const CerenkovHistogramFile & histFile, int numberOfPoints)
{
  return CerenkovAngularDistribution::GetLemoineHistogram(histFile.GetShowerZenith(), numberOfPoints);
}

TH1 * CerenkovAngularDistribution::GetGillerHistogram( const CerenkovHistogramFile & histFile, int numberOfPoints)
{
  return CerenkovAngularDistribution::GetGillerHistogram(histFile.GetMeanHeight(),histFile.GetAge(), numberOfPoints);
}

TH1 * CerenkovAngularDistribution::GetNerlingHistogram(const CerenkovHistogramFile & histFile, int numberOfPoints)
{
  return CerenkovAngularDistribution::GetNerlingHistogram(histFile.GetAge(),0.36133/std::sqrt(histFile.GetRefractiveIndex()-1), numberOfPoints);
}
