#pragma once

#include <vector>

#include <TH1.h>
#include <Fit/FitResult.h>

#include <cerenkovHistogramFile.h>

/**
  @brief A class to compute the angular distribution function given parameters.

  This class uses the proposed function to compute the angular distribution function. Parameters can be set manually or
  be estimated from a given histogram by using a CerenkovFitter object for a single histogram.

  @author Luan Arbeletche
  @date July 2020
 */
class CerenkovAngularDistribution
{
private:
  double emissionAngle; /**< Value of theta_em (maximum emission angle of Cherenkov photons) */
  double parTheta_1; /**< Value of the parameter theta_1 */
  double parRatio; /**< Value of the parameter theta_2/theta_1 */
  double parNu; /**< Value of the parameter Nu */
  double parEps; /**< Value of the parameter Epsilon */
  double errTheta_1; /**< Uncertainty on the parameter Theta_1 */
  double errRatio; /**< Uncertainty on the parameter Theta_2/Theta_1*/
  double errNu; /**< Uncertainty on the parameter Nu */
  double errEps; /**< Uncertainty on the parameter Epsilon */
  ROOT::Fit::FitResult fitResult; /**< The fit result */
  
  std::vector<double> parameterSeedArray; /**< Parameter seeds to perform the fit */
  std::vector<double> lowerParameterBound; /**< Array of lower parameter bounds */
  std::vector<double> upperParameterBound; /**< Array of upper parameter bounds */
  
  /**
    Function used to compute the parameters given shower age, refractive index - 1 and shower energy. Used when creating a CerenkovFitter object
    within the EstimateParameters() function.
    @param functionParameters - Output array with the parameters used to compute the angular distribution function
    @param p - Input array with the fit coefficients
    @param age - Shower age
    @param delta - Refractive index minus one
    @param energy - Shower energy
   */ 
  static void ParameterCalculator(double * functionParameters, const double *& p, const double &age, const double &delta, const double &energy);

public:
  /**
    Constructor.
   */ 
  CerenkovAngularDistribution();

  /**
    Estimate parameters of the angular distribution function to fit the given histogram.
    @param cerHistogramFileRef - CherenkovHistogramFile object from which the data histogram will be extracted
    @return The ROOT::Fit::FitResult object obtained in the fit.
   */ 
  ROOT::Fit::FitResult EstimateParameters(const CerenkovHistogramFile & cerHistogramFileRef);

  /**
    Manually set the parameters of the angular distribution function.
    @param tch - The maximum cherenkov emission angle (equal to theta_em)
    @param t1 - Value of THeta_1
    @param r - value of Theta_2 / Theta_1
    @param nu - Value of Nu
    @param eps - Value of Epsilon
   */ 
  void SetParameters(double tch, double t1, double r, double nu, double eps);

  /**
    Simple getter to retrieve the fit result if EstimateParameters() has been called
    @return The ROOT::Fit::FitResult object.
   */ 
  ROOT::Fit::FitResult GetFitResult() const {
    return this->fitResult;
  }

  /**
    Evaluate the angular distribution function at a given number of points within the specified theta interval. Optionally, scales
    the angular distribution function by the given value.
    @param numberOfPoints - Number of points to compute the angular distribution function
    @param minTheta - Low edge of theta interval
    @param maxTheta - Upper edge of theta interval
    @param normalization - (optional) Scale factor to apply to the computed function
    @return A vector containing the computed values in ascending order of theta
   */ 
  std::vector<double> EvalVector(int numberOfPoints, double minTheta, double maxTheta, double normalization = 1) const;

  /**
    Similar to EvalVector(), but returns a histogram instead of a vector with function values.
    @param numberOfPoints - Number of points to compute the angular distribution function
    @param minTheta - Low edge of theta interval
    @param maxTheta - Upper edge of theta interval
    @param normalization - (optional) Scale factor to apply to the computed function
    @return A histogram with the computed function values
   */ 
  TH1 * GetFitHistogram(int numberOfPoints, double minTheta, double maxTheta, double normalization = 1) const;
  
  /**
    Function to set the seeds used in the fit procedure.
    @param t1 - Initial value of Theta_1
    @param r - Initial value of Theta_2 / Theta_1
    @param nu - Initial value of Nu
    @param eps - Initial value of Epsilon
   */ 
  void SetParameterSeed(double t1, double r, double nu, double eps);

  /**
    Function to set the lower bounds used in the fit procedure.
    @param t1 - Initial value of Theta_1
    @param r - Initial value of Theta_2 / Theta_1
    @param nu - Initial value of Nu
    @param eps - Initial value of Epsilon
   */ 
  void SetLowerBounds(double t1, double r, double nu, double eps);

  /**
    Function to set the upper bounds used in the fit procedure.
    @param t1 - Initial value of Theta_1
    @param r - Initial value of Theta_2 / Theta_1
    @param nu - Initial value of Nu
    @param eps - Initial value of Epsilon
   */ 
  void SetUpperBounds(double t1, double r, double nu, double eps);

  //
  // Below are static methods to compute the distributions from other references
  //

  /**
    Computes the angular distribution function as published by Lemoine-Goumard et al.
    @param showerZenith - Shower zenith angle
    @param numberOfPoints - (optional) (default 600) Number of points in the output histogram
    @return A histogram with the angular distribution
   */ 
  static TH1 * GetLemoineHistogram(double showerZenith, int numberOfPoints = 600);

  /**
    Computes the angular distribution function as published by Giller et al.
    @param height - Altitude
    @param age - Shower age
    @param numberOfPoints - (optional) (default 600) Number of points in the output histogram
    @return A histogram with the angular distribution
   */ 
  static TH1 * GetGillerHistogram(double height, double age, int numberOfPoints = 600);

  /**
    Computes the angular distribution function as published by Nerling et al.
    @param age - Shower age
    @param threshold - Energy threshold for the emission of Cherenkov light
    @param numberOfPoints - (optional) (default 600) Number of points in the output histogram
    @return A histogram with the angular distribution
   */ 
  static TH1 * GetNerlingHistogram(double age, double threshold, int numberOfPoints = 600);

  /**
    Computes the angular distribution function as published by Lemoine-Goumard et al.
    @param histFile - A CerenkovHistogramFile object from which shower information will be extracted
    @return A histogram with the angular distribution
   */ 
  static TH1 * GetLemoineHistogram(const CerenkovHistogramFile & histFile, int numberOfPoints = 600);

  /**
    Computes the angular distribution function as published by Giller et al.
    @param histFile - A CerenkovHistogramFile object from which shower information will be extracted
    @return A histogram with the angular distribution
   */ 
  static TH1 * GetGillerHistogram( const CerenkovHistogramFile & histFile, int numberOfPoints = 600);

  /**
    Computes the angular distribution function as published by Nerling et al.
    @param histFile - A CerenkovHistogramFile object from which shower information will be extracted
    @return A histogram with the angular distribution
   */ 
  static TH1 * GetNerlingHistogram(const CerenkovHistogramFile & histFile, int numberOfPoints = 600);
};
