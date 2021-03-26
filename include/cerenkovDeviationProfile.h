#pragma once

// STL dependencies
#include <string>
#include <map>
// ROOT dependencies
#include <TProfile.h>
// Internal dependencies
#include <cerenkovHistogramFile.h>
#include <cerenkovAngularDistribution.h>

/**
  @brief A class to compute profiles with the average relative deviation between a parametrization and the MC data.

  This class is intended to encapsulate the computation of the relative deviation between the proposed parametrization and the MC
  data used to obtain the parametrization. It provides methods to write tables with computed data to be used by an external program
  to plot the deviation profiles.

  @author Luan Arbeletche
  @date July 2020
 */
class CerenkovDeviationProfile
{
private:
  std::string name; /**< Name of this deviation profile */
  int numberOfPoints; /**< Number of points in the profile */
  double minTheta; /**< Min value for x axis*/
  double maxTheta; /**< Max value for x axis*/
  double maxThetaPlot; /**< Max value of x axis when printing data*/
  
  std::map<int, TProfile> mapOfProfiles; /**< A structure to hold the deviation profiles. The index is int(showerAge * 10) */
  
  /**
    Checks if the shower age associated with the given histogram matches some of the profiles in the profile map mapOfProfiles.
    @param histogramFile - CerenkovHistogramFile object from which the histogram will be consulted
    @return Index of mapOfProfiles corresonding with the age of this histogram, if found. Otherwise, return -1.
   */
  int CheckAge(const CerenkovHistogramFile & histogramFile);

public:
  /**
    Constructor.
    @param deviationProfileName - Name of the deviation profile
    @param npt - Number of points in the profiles
    @param min - Minimum value of the angle theta
    @param max - Maximum value of the angle theta
   */
  CerenkovDeviationProfile(std::string deviationProfileName, int npt, double min, double max);

  /**
    Destructor.
   */
  ~CerenkovDeviationProfile();

  /**
    Creates a new, empty profile for the given value of shower age.
    @param age - Shower age associated with the new deviation profile
   */
  void AddNewProfile(double age);

  /**
    Adds deviation data from the given combination of data histogram and angular distribution function.
    @param histogramFile - Histogram file from which the reference data (MC) will be consulted.
    @param angularDist - Angular distribution describing the histogram data.
   */
  void AddData(const CerenkovHistogramFile & histogramFile, const CerenkovAngularDistribution & angularDist);

  /**
    Adds deviation data from the given combination of data histogram and angular distribution function given as a histogram.
    @param histogramFile - Histogram file from which the reference data (MC) will be consulted.
    @param inputFitHistogram - Histogram with the computed angular distribution.
   */
  void AddData(const CerenkovHistogramFile & histogramFile, TH1 * inputFitHistogram);

  /**
    Plots and prints all computed profiles to a .pdf file in the given path.
    @param outputPath - Path to the output .pdf file with plots.
   */
  void Print(std::string outputPath);

  /**
    Write a table with data from the computed profiles in the given output file.
    @param outputPath - Path to the output file
    @param rebinFactor - (optional) Rebin deviation profiles if != 1
   */
  void WriteTable(std::string outputPath, int rebinFactor = 1);

  /**
    Sets a maximum value of theta on the printed deviation profiles.
    @param value - New value for the maximum theta.
   */
  void SetMaxThetaPlot(double value){
    this->maxThetaPlot = value;
  }
};

