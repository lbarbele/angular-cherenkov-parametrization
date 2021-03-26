#pragma once

// STL dependencies
#include <map>
#include <string>
#include <iostream>
// ROOT dependencies
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
// Internal dependencies
#include <cerenkovHistogramFile.h>
#include <cerenkovAngularDistribution.h>

/**
  @brief A class to standardize plots generated in this study.

  This class derives from TCanvas and aims to provide standard plots of the angular distribution of Cherenkov photons, both from simulations
  and the parametrized ones, and also to build data tables that can be used by an external plot program (like pgfplots).

  @author Luan Arbeletche
  @date July 2020
 */
class CerenkovCanvas : public TCanvas
{
private:
  std::string outputFileName; /**< Name of the output file */
  
  double minTheta; /**< Left edge of x axis */
  double midTheta; /**< Threshold separating large amd small angular regions (to build inset plots) */
  double maxTheta; /**< Upper edge of x axis */
  
  std::unique_ptr<TH2> hAxisLin; /**< Pointer to a histogram defining axes for the linear plot in the small theta region */
  std::unique_ptr<TH2> hAxisLog; /**< Pointer to a histogram defining axes for the log scale plot in large theta region */
  std::string drawOptionHistogram; /**< Options to be used in the Draw() method of the MC histograms */
  std::string drawOptionFunction; /**< Options to be used in the Draw() method of the parametrized/fitted distributions */
  
  std::map<std::string, std::map< std::string, TH1D>  > mapOfData; /**< A map of data (used to build tables) storing values according to keys [tableName][columnName] */

public:
  /**
    Main constructor.
    Creates a canvas object with a name for the output file
    @param outputFileName - path to the output .pdf file
   */
  CerenkovCanvas(std::string outputFileName);

  /**
    Destructor.
   */
  ~CerenkovCanvas();
  
  /**
    Adds a data histogram and an angular distribution function and send the page to the output .pdf file.
    @param histFilePtr - Pointer to a CerenkovHistogramFile object from which the histogram will be extracted
    @param angularDistPtr - Pointer to a CerenkovAngularDistribution object to be drawn
    @param doNormalize - (optional) Normalize histogram and function, if activated
   */
  void AddPage(const CerenkovHistogramFile * histFilePtr, const CerenkovAngularDistribution * angularDistPtr, bool doNormalize = true); 
  
  /**
    Add a generic TH1* to the current page and optionally print the page
    @param histogramPtr - The histogram to be plotted in the current page
    @param doOpenPage - (optional) Flag telling if a new page is to be created
    @param doClosePage - (optional) Flag telling if page is to be closed after ploting this histogram
    @param doNormalize - (optional) Flag asking for normalization
    @param label - (optional) Label of the current page
   */
  void AddData(TH1 * histogramPtr, bool doOpenPage = true, bool doClosePage = true, bool doNormalize = false, std::string label = "");
  
  // Add histogram as column "columnName" to table called "tableName"

  /**
    Add data from the histogram passed as input to the table tableName as a column named columName
    @param tableName - Name of the table to which data is to be added
    @param columnName - Name of the column in this table
    @param histFilePtr - Pointer to a CerenkovHistogramFile from which a histogram will be read
   */
  void MakeDataTable(std::string tableName, std::string columnName, const CerenkovHistogramFile * histFilePtr);

  /**
    Add data from the angular distribution function passed as input to the table tableName as a column named columName
    @param tableName - Name of the table to which data is to be added
    @param columnName - Name of the column in this table
    @param angularDistPtr - Pointer to a CerenkovAngularDistribution object from which data will be computed
    @param numberOfPoints - Number of points to compute and add to the table
   */
  void MakeDataTable(std::string tableName, std::string columnName, const CerenkovAngularDistribution * angularDistPtr, int numberOfPoints = 600);

  /**
    Add data from the histogram passed as input to the table tableName as a column named columName
    @param tableName - Name of the table to which data is to be added
    @param columnName - Name of the column in this table
    @param histFilePtr - Histogram from which data will be extracted
   */
  void MakeDataTable(std::string tableName, std::string columnName, TH1 * inputHistogram);
                    
  /**
    Print a table whose name is passed as input to a given stream, optionally setting a cut on the maximum value of theta.
    @param tableName - Name of the internal table to be printed
    @param outStream - Output stream to which the data will be sent as a formatted table
    @param maxPrintThetaDeg - (optional) Sets a maximum value of theta on the output data
   */ 
  void PrintTable(std::string tableName, std::ostream && outStream, double maxPrintThetaDeg = 60);
  
  /**
    Sets options to be sent to the ROOT Draw() function when drawing histograms.
    @param sopt - The new drawing options list
   */
  void SetHistDrawOption(std::string sopt) 
  {
    this->drawOptionHistogram = sopt + " same";
  }
  
  /**
    Sets options to be sent to the ROOT Draw() function when drawing functions.
    @param sopt - The new drawing options list
   */
  void SetFuncDrawOption(std::string sopt)
  {
    this->drawOptionFunction = sopt + " same";
  }
};
