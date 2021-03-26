#include <cmath>
#include <algorithm>
#include <memory>
#include <iomanip>
#include <iostream>

#include <TH2.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TSystem.h>

#include <cerenkovCanvas.h>

CerenkovCanvas::CerenkovCanvas(std::string outputFileNameInput) :
  TCanvas("cerCanvas","",1200,600),
  outputFileName(outputFileNameInput),
  minTheta(0),
  midTheta(5.0*std::acos(-1)/180.0),
  maxTheta(60.0*std::acos(-1)/180.0),
  drawOptionHistogram("hist same"),
  drawOptionFunction("hist l same"),
  hAxisLin(nullptr),
  hAxisLog(nullptr)
{
  // Divide the canvas
  this->Divide(2,1);
  
  // Open the output file
  this->Print((this->outputFileName + "[").c_str());
}

CerenkovCanvas::~CerenkovCanvas()
{
  // Close the output file
  this->Print((this->outputFileName + "]").c_str());
}

void CerenkovCanvas::AddPage(const CerenkovHistogramFile * histFilePtr, const CerenkovAngularDistribution * angularDistPtr, bool doNormalize)
{
  // If histogram is not present, exit
  if (!histFilePtr) {
    return;
  }
    
  // Get pointer to histogram
  TH1 * hist = histFilePtr->GetThetaHistogramPointer();

  if (doNormalize) {
    hist->Scale(1.0/hist->Integral("width"));
  }

  // Build a graph of the fit function, if angularDistPtr is given
  std::unique_ptr<TH1> hFit( nullptr );

  if (angularDistPtr) {
    hFit.reset(angularDistPtr->GetFitHistogram(hist->GetNbinsX(),this->minTheta,this->maxTheta, hist->Integral()));
  }
    
  // Parametrization by Lemoine-Goumard et al
  std::unique_ptr<TH1> hLemoine( CerenkovAngularDistribution::GetLemoineHistogram(*histFilePtr,hist->GetNbinsX()) );
  
  // Parametrization by Giller et al
  std::unique_ptr<TH1> hGiller( CerenkovAngularDistribution::GetGillerHistogram(*histFilePtr,hist->GetNbinsX()) );
  
  // Parametrization by Nerling et al
  std::unique_ptr<TH1> hNerling( CerenkovAngularDistribution::GetNerlingHistogram(*histFilePtr,hist->GetNbinsX()) );

  // Build axes
  double yMax = hist->GetMaximum();
  double yMin = hist->GetMinimum(0);
  double delta = std::exp(0.05*(std::log(yMax) - std::log(yMin)));
  
  this->hAxisLin.reset( new TH2D("haxislin","",100,this->minTheta,this->midTheta,100,0.0,yMax*1.05) );
  this->hAxisLog.reset( new TH2D("haxislog","",100,this->minTheta,this->maxTheta,100,yMin/delta,yMax*delta) );
  
  this->hAxisLog->SetStats(false);
  this->hAxisLin->SetStats(false);
  
  // Histogram colors
  hist->SetLineColor(kBlack);

  if (hFit) {
    hFit->SetLineColor(kViolet-3);
  }

  if (hLemoine) {
    hLemoine->SetLineColor(kGreen+4);
  }

  if (hGiller) {
    hGiller->SetLineColor(kRed);
  }
  
  if (hNerling) {
    hNerling->SetLineColor(kBlue);
  }
  
  // Draw left page
  this->cd(1);
  gPad->Clear();
  
  this->hAxisLin->Draw("axis");
  hist->Draw(this->drawOptionHistogram.c_str());

  if (hFit) {
    hFit->Draw(this->drawOptionFunction.c_str());
  }

  if (hLemoine) {
    hLemoine->Draw(this->drawOptionFunction.c_str());
  }
  
  if (hGiller) {
    hGiller->Draw(this->drawOptionFunction.c_str());
  }

  if (hNerling) {
    hNerling->Draw(this->drawOptionFunction.c_str());
  }
  
  TLatex labelWithShowerAge;
  std::string ageString = std::to_string(histFilePtr->GetAge()).substr(0,4);
  labelWithShowerAge.DrawLatexNDC(0.05,0.95,("s = " + ageString).c_str());
  
  // Draw right page
  this->cd(2);
  gPad->Clear();
  gPad->SetLogy();
  
  this->hAxisLog->Draw("axis");
  hist->Draw(this->drawOptionHistogram.c_str());

  if (hFit) {
    hFit->Draw(this->drawOptionFunction.c_str());
  }

  if (hLemoine) {
    hLemoine->Draw(this->drawOptionFunction.c_str());
  }

  if (hGiller) {
    hGiller->Draw(this->drawOptionFunction.c_str());
  }

  if (hNerling) {
    hNerling->Draw(this->drawOptionFunction.c_str());
  }
  
  // Print this page
  this->Print(this->outputFileName.c_str());
  
  // Reset the axis histogram pointers
  this->hAxisLin.reset(nullptr);
  this->hAxisLog.reset(nullptr);
}

void CerenkovCanvas::AddData(TH1 * histogramPtr, bool doOpenPage, bool doClosePage, bool doNormalize, std::string label)
{
  // Check if *histogramPtr exists
  if (!histogramPtr) {
    return;
  }
    
  // Check normalizations
  if (doNormalize) {
    histogramPtr->Scale(1.0/histogramPtr->Integral("width"));
  }
    
  // If this is a new page, create the axes, clear the canvas, draw the axes
  if (doOpenPage) {
   
    double yMax = histogramPtr->GetMaximum();
    double yMin = histogramPtr->GetMinimum(0);
    double delta = std::exp(0.05*(std::log(yMax) - std::log(yMin)));
  
    this->hAxisLin.reset( new TH2D("haxislin","",100,this->minTheta,this->midTheta,100,0.0,yMax*1.05) );
    this->hAxisLog.reset( new TH2D("haxislog","",100,this->minTheta,this->maxTheta,100,yMin/delta,yMax*delta) );
    
    this->hAxisLog->SetStats(false);
    this->hAxisLin->SetStats(false);
    
    // Clear and draw axis on left page
    this->cd(1);
    gPad->Clear();
    this->hAxisLin->Draw("axis");
    
    // If there is a label, draw it
    if (label.size() > 0) {
      TLatex labelWithShowerAge;
      labelWithShowerAge.DrawLatexNDC(0.05,0.95,label.c_str());
    }
    
    // Do the same on the right page
    this->cd(2);
    gPad->Clear();
    this->hAxisLog->Draw("axis");
    gPad->SetLogy();
  }
  
  // Draw histogram, if pointer exists
  if (histogramPtr) {
    this->cd(1);
    histogramPtr->Draw(this->drawOptionHistogram.c_str());
  
    this->cd(2);
    histogramPtr->Draw(this->drawOptionHistogram.c_str());
  }
  
  // If page is to be closed, do it
  if (doClosePage) {
    // Print page
    this->Print(this->outputFileName.c_str());
    
    // Reset axes
    this->hAxisLin.reset(nullptr);
    this->hAxisLog.reset(nullptr);
  }
}

void CerenkovCanvas::MakeDataTable(std::string tableName, std::string columnName, const CerenkovHistogramFile * histFilePtr)
{
  // Check if input pointer is valid
  if (!histFilePtr) {
    return;
  }
    
  // Check if column name is available
  if (this->mapOfData[tableName].find(columnName) != this->mapOfData[tableName].end()) {
    std::cerr << "CerenkovCanvas::MakeDataTable: a column called \"" << columnName << "\" already exists in table \"" << tableName << "\"." << std::endl;
    return;
  }
  
  // Put data on the map
  this->mapOfData[tableName].emplace(columnName, *(TH1D*)histFilePtr->GetThetaHistogramPointer() );
}

void CerenkovCanvas::MakeDataTable(std::string tableName, std::string columnName, const CerenkovAngularDistribution * angularDistPtr, int numberOfPoints)
{
  // Check if input pointer is valid
  if (!angularDistPtr) {
    return;
  }
    
  // Check if column name is available
  if (this->mapOfData[tableName].find(columnName) != this->mapOfData[tableName].end()) {
    std::cerr << "CerenkovCanvas::MakeDataTable: a column called \"" << columnName << "\" already exists in table \"" << tableName << "\"." << std::endl;
    return;
  }
  
  // Put data on map
  std::unique_ptr<TH1D> hist( (TH1D*) angularDistPtr->GetFitHistogram(numberOfPoints, 0, std::acos(-1)/3.0, std::acos(-1)/(3.0*numberOfPoints)) );
  this->mapOfData[tableName].emplace(columnName, *hist);
}

void CerenkovCanvas::MakeDataTable(std::string tableName, std::string columnName, TH1 * inputHistogram)
{
  // Check if input pointer is valid
  if (!inputHistogram) {
    return;
  }
    
  // Check if column name is available
  if (this->mapOfData[tableName].find(columnName) != this->mapOfData[tableName].end()) {
    std::cerr << "CerenkovCanvas::MakeDataTable: a column called \"" << columnName << "\" already exists in table \"" << tableName << "\"." << std::endl;
    return;
  }
  
  // Put data on map
  this->mapOfData[tableName].emplace(columnName, *(TH1D*)inputHistogram);
}

void CerenkovCanvas::PrintTable(std::string tableName, std::ostream && outStream, double maxPrintThetaDeg)
{
  // normalize data (this could be optional, but for now it is not)
  for (auto & ppp : this->mapOfData[tableName]) {
    ppp.second.Scale(1.0/ppp.second.Integral("width"));
  }
  
  // print table header (column names)
  outStream << std::setw(13) << "theta";
  outStream << std::setw(13) << "thetaLeft";

  for (auto & ppp : this->mapOfData[tableName]) {
    outStream << std::setw(13) << ppp.first;
  }
  outStream << std::endl;
  
  // get a reference histogram used to print the abscissas and get number of bins
  auto & referenceHistogram = (*this->mapOfData[tableName].begin()).second;
  
  // loop over bins and send to output stream
  for (int iBin = 1; iBin <= referenceHistogram.GetNbinsX(); iBin++) {

    double theta = referenceHistogram.GetBinCenter(iBin)*57.295779513;
    
    if (theta > maxPrintThetaDeg)
      break;
    
    outStream << std::setw(13) << theta;
    outStream << std::setw(13) << referenceHistogram.GetBinLowEdge(iBin)*57.295779513;
    for (auto & ppp : this->mapOfData[tableName]) outStream << std::setw(13) << ppp.second.GetBinContent(iBin);
    outStream << std::endl;
  }
}
