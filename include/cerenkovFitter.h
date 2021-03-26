#pragma once

// STL dependencies
#include <vector>
#include <mutex>
// ROOT dependencies
#include <TH1.h>
#include <Fit/FitResult.h>
// Internal dependencies
#include <cerenkovHistogramFile.h>

/**
  @brief Class to parametrize the angular distribution of Cherenkov photons.

  This class is aimed to build a parametrization of the angular distribution of Cherenkov photons using the angular distribution function
  proposed in 10.1140/epjc/s10052-021-08971-7, but allowing to vary the functional form in which the function parameters evolve with
  shower age, refractive index, and primary energy. This functional form is passed as a function via SetParameterFunction(). Besides encapsulating
  the objective function used in the minimization procedure, this function encapsulates data in the form of histograms. These histograms are
  added to the fit procedure via calls to AddHistogram().

  An usage example can be found in fit.cpp.

  @author Luan Arbeletche
  @date July 2020
 */
class CerenkovFitter
{
private: 

  /**
    @brief A class to hold a histogram of the angular distribution of Cherenkov photons and other information used to compute the deviance function.
    @author Luan Arbeletche
    @date July 2020
   */
  struct DataEntry
  { 
    std::vector<double> vector_y; /**< Auxiliar quantity. See constructor for definition.*/
    std::vector<double> vector_x; /**< Auxiliar quantity. See constructor for definition. */
    std::vector<double> vector_j; /**< Auxiliar quantity. See constructor for definition. */
    double A; /**< Auxiliar quantity. See constructor for definition.*/
    double B; /**< Auxiliar quantity. See constructor for definition.*/
    double C; /**< Auxiliar quantity. See constructor for definition.*/
    double N; /**< Auxiliar quantity. See constructor for definition.*/
    double showerEnergy; /**< Shower energy */
    double showerZenith; /**< Shower zenith */
    double refractiveIndex; /**< Refractive index*/
    double showerAge; /**< Shower age*/
    double emissionAngle; /**< Maximum emission angle of cherenkov light*/
    double minTheta; /**< Lower edge of histogram abscissa */
    double maxTheta;/**< Upper edge of histogram abscissa */
    int numberOfBins;/**< Number of bins in the histogram */

    DataEntry * previousEntry; /**< Previous entry in the list */
    DataEntry * nextEntry; /**< Next entry in the list */
    
    /**
      Constructor.
      @param energyInput - Shower energy
      @param zenithInput - Zenith angle
      @param ageInput - Shower age
      @param refIndexInput - Refractive index
      @param histogram - Reference to the histogram that this entry will hold
      @param previousInput - Pointer to a DataEntry object after which this entry will be placed
     */
    DataEntry(double energyInput, double zenithInput, double ageInput, double refIndexInput, const TH1 & histogram, DataEntry * previousInput = nullptr);

    /**
      Constructor.
     */
    ~DataEntry();
  };
  
  int numberOfParameters; /**< Number of free parameters in the parametrization */
  double minAge; /**< Minimum age of data histograms to be considered */
  double maxAge; /**< Maximum age of data histograms to be considered */
  bool isVerbose; /**< Verbosity flag */

  DataEntry * headOfEntryList; /**< First histogram on the data list */
  DataEntry* tailOfEntryList; /**< Last histogram on the data list */
  void (*parCalcPtr)(double *,
                     const double *&,
                     const double&,
                     const double&,
                     const double&); /**< Pointer to a function given by the user that is used to compute parameters in the minimization procedure */
  
  DataEntry * nextEntryFromList; /**< A pointer to the next entry, must be read and modified from GetNextEntryFromListAsync() */
  std::mutex nextEntryReaderBarrier; /**< barrier to avoid that multi threads access/modify nextEntryFromList at the same time */
  
  /**
    Compute deviance function asynchronously using a given array of parameters. May be called multiple times from different threads.
    @param parameterArray - Array of parameters used to compute the objective function
    @return Value of resulting deviance (likelihood function)
   */
  double ComputeDevianceAsync(const double * parameterArray);

  /**
    Asynchronously get the pointer to the next data entry and modify the internal pointer.
    @return Pointer to next data entry.
   */
  DataEntry * GetNextEntryFromListAsync(); // Return nextEntryFromList and modify it 
  
public:
  /**
    Default constructor. Needs no parameters.
   */
  CerenkovFitter();

  /**
    Adds a histogram to the list of entries optionally normalizing it. If given normalization factor is present and is != 1, scale the histogram by this value.
    @param cerHistogramFileRef - File from which the histogram will be read using GetThetaHistogramPointer()
    @param normalizationFactor - Factor to which the histogram will be scaled
    @return True if histogram is within specified age limits and was added, otherwise false.
   */
  bool AddHistogram(const CerenkovHistogramFile & cerHistogramFileRef, double normalizationFactor = 1);

  /**
    Function to actually perform the parametrization.
    @param seedVector - Initial values to the fit parameters
    @param lowerVector - (optional) Lower bounds to fit parameters
    @param upperVector - (optional) Upper bounds to fit parameters
    @return A FitResult object containing information from the minimization procedure.
   */
  ROOT::Fit::FitResult Fit(std::vector<double> seedVector, std::vector<double> lowerVector = std::vector<double>(0), std::vector<double> upperVector = std::vector<double>(0));

  /**
    This overloaded operator () is responsible to compute the deviance function.
    @param seedVector - Initial values to the fit parameters
    @param lowerVector - (optional) Lower bounds to fit parameters
    @param upperVector - (optional) Upper bounds to fit parameters
    @return A FitResult object containing information from the minimization procedure.
   */
  double operator()(const double * parameters);
  
  /**
    Function to set the verbosity flag.
    @param wantVerbosity - New value of the verbosity flag.
   */
  void SetVerbosity(bool wantVerbosity)
  {
    this->isVerbose = wantVerbosity;
  }
  
  /**
    Function to set the function that will be used in the parametrization to compute the dependent parameters in terms of the free parameters.
    @param inputFunction - A function with signature void(double *, const double *&, const double&, const double&, const double&)
   */
  void SetParameterFunction(void (*inputFunction)(double *, const double *&, const double&, const double&, const double&))
  {
    this->parCalcPtr = inputFunction;
  }
  
  /**
    Function to set the age interval within which histograms will be considered in the deviance function.
    @param minAgeInput - New value for the minimum age
    @param maxAgeInput - New value for the maximum age
   */
  void SetAgeInterval(double minAgeInput, double maxAgeInput)
  {
    this->minAge = minAgeInput;
    this->maxAge = maxAgeInput;
  }
};
