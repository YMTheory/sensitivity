#ifndef nEXOUtil_hh
#define nEXOUtil_hh

#include <iostream>
#include <vector>

#include "TH2.h"
#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TList.h"
#include "TH2D.h"
#include "TF1.h"
#include "Math/BrentRootFinder.h"
#include "Math/WrappedTF1.h"
#include "Math/ProbFuncMathCore.h"
#include "TGaxis.h"
#include "TMath.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TMatrixDSym.h"
#include "TTree.h"
#include "TLegend.h"
#include "TChain.h"

#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooAbsReal.h"
#include "RooGenericPdf.h"
#include "RooGlobalFunc.h"
#include "RooFitResult.h"


enum MATRIX_ELEMENT {GCM=0, QRPA2, RQRPA, NSM, IBM2};

class nEXOUtil : public TObject
{
public:
  nEXOUtil();
  virtual ~nEXOUtil();
  
private:
  Double_t 	fMatElem[5];
  TString	fMatElemName[5];
  TString	fMatElemTitle[5];
  TString	fMatElemAxisTitle[5];
  TString	fXaxisTitle;
  TString	fYaxisTitle;
  Double_t 	fDiscPlotRangeX[2];
  Double_t 	fDiscPlotRangeY[2];
  static constexpr Int_t	fNBinsX = 2800;
  static constexpr Int_t	fNBinsY = 650;
  Double_t	fFidVolXY;
  Double_t	fFidVolZZ;
  static constexpr Double_t	fXeMass = 4779.96;
  static constexpr Double_t	fEnrichment = 0.9;
  static constexpr Double_t	fXeMolarMass = 0.136;
  static constexpr Double_t	fAvogadrosNumber = 6.022e23;
  
  Double_t getMatElem(enum MATRIX_ELEMENT matElem) 			{return fMatElem[matElem];}
  TString getMatElemName(enum MATRIX_ELEMENT matElem) 		{return fMatElemName[matElem];}
  TString getMatElemTitle(enum MATRIX_ELEMENT matElem) 		{return fMatElemTitle[matElem];}
  TString getMatElemAxisTitle(enum MATRIX_ELEMENT matElem) {return fMatElemAxisTitle[matElem];}

public:
  //median counts calculation
  void PlotBb0nUL(TString inFileName, TString outFileName, Int_t withBb0n=false);
  double getMedian(TH1D* h1, double& lo, double& hi);
  double getHalfLife(TString inFileName, Double_t yrs, Double_t xeMass);
  double getBkgdCounts(TString inFileName, int vol);
  
  //histogram utility functions
  Double_t getFractional2DIntegral(TH2* h, Double_t xLo, Double_t xHi, Double_t yLo, Double_t yHi);
  TH2D* combineHistos(std::vector<TH2D*>* histos, std::vector<Double_t>* fracs, Double_t norm=1.);
  void rebin2DHistoFile(const char* inFileName, const char* outFileName, Int_t rebinX, Int_t rebinY);
  void normalizeHistoFile(const char* inFileName, const char* outFileName);
  void smoothHisto(TH2* histo, Bool_t isSS, Int_t startBin);
  void getSmoothArray(Double_t* smoothArray, TH2* histo, Bool_t isSS);
  Double_t getNumXeAtoms(Double_t xeMass) {return xeMass*fEnrichment*fAvogadrosNumber/fXeMolarMass;}

  //discovery potential calculating functions
  Double_t getLeastDetectableSignal(Double_t prob, Double_t nBkgd, Double_t nCrit);
  Double_t getCriticalCounts(Double_t alpha, Double_t nBkgd);
  Double_t getDiscoveryCounts(Double_t bkgCounts, Double_t gausSignif, Double_t prob, Double_t eff=1.);
  std::vector<Double_t> getListOfDiscoveryCounts(std::vector<Double_t> bkgCounts, Double_t gausSignif, Double_t prob, Double_t eff);
  Double_t getBkgdCountsInRange(const char* histFile, Double_t yrsExp, Int_t nPdfs, Double_t meanPerYear[], TString pdfNames[], Double_t* xRange, Double_t* yRange, Bool_t isBaTag=false);
  Double_t getHalfLifeForCounts(Double_t nYrs, Double_t counts, Double_t xeMass);
  Double_t getCutValue(TH1D* histo, Double_t gausSignif);
  
  //drawing functions
  TGraph* getInvertedGraph(enum MATRIX_ELEMENT matElem);
  TGraph* getNormalGraph(enum MATRIX_ELEMENT matElem);
  TCanvas* getEmptyDiscoveryPlot(enum MATRIX_ELEMENT matElem, Double_t xMin=0., Double_t xMax=10., Double_t yMin=1.e25, Double_t yMax=2.e29);
  
  //other functions
  Double_t getCo60CountsAfterYears(Double_t meanPerYear, Double_t yrs); 
  
  ClassDef(nEXOUtil,1)	// Class with utility functions for nEXO sensitivity
};

#endif
