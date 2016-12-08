#ifndef nEXOUtils_hh
#define nEXOUtils_hh

#include <iostream>
#include <vector>

#include "TString.h"
#include "TChain.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TFile.h"
#include "Math/BrentRootFinder.h"
#include "Math/WrappedTF1.h"
#include "Math/ProbFuncMathCore.h"
#include "TFeldmanCousins.h"

#include "nEXONuclearMatrixElement.hh"
#include "nEXONuOscPars.hh"

namespace nEXOUtils
{
  // Utility functions and types for the nEXO sensitivity
  
  // useful constants
  static const Double_t	XE136_MOLAR_MASS = 0.136;
  static const Double_t	AVOGADRO_NUMBER = 6.022e23;

  // functions for sensitivity calculations
  Double_t GetNumXeAtoms(Double_t xeMass, Double_t enrichment = 0.9);
  Double_t GetHalfLife(Double_t counts, Double_t yrs, Double_t xeMass);
  Double_t GetSensHalfLife(TString inFileName, Double_t yrs, Double_t xeMass, TString signalName = "signal", double eff = 1.);
  
  Double_t GetBackgroundCounts(TString inFileName, TString varName = "bkg_tot", TString treeName = "tree", bool useMean = false);  
  Double_t GetBackgroundCounts(Long64_t length, Double_t* counts, bool useMean = false);
  Double_t GetDiscoveryCutValue(TH1D* histo, Double_t gausSignif);
  Double_t GetDiscoveryCutValue(TH1D* histo, Double_t gausSignif);
  
  Double_t GetDiscHalfLife(Double_t prob, Double_t sigma, Long64_t length, Double_t* counts, Double_t yrs, const char* filename, double mass, double max = 450., double min = -50., int bins = 10000, const char* treename = "tree", const char* cut = "stat_sig==0 && stat_bkg==0 && covQual_sig==3 && covQual_bkg==3");
  Double_t GetDiscoveryCounts(Double_t prob, Double_t sigma, Long64_t length, Double_t* counts, Double_t yrs, const char* filename, double max = 450., double min = -50., int bins = 10000, const char* treename = "tree", const char* cut = "stat_sig==0 && stat_bkg==0 && covQual_sig==3 && covQual_bkg==3"); //"fitres_sig._status==0 && fitres_bkg._status==0 && fitres_sig._covQual==3 && fitres_bkg._covQual==3"
  bool BuildNLLRatioHistogram(TH1D& histo, const char* filename, const char* treename = "tree", const char* cut = "stat_sig==0 && stat_bkg==0 && covQual_sig==3 && covQual_bkg==3");

  // Counting experiment
  // Sensitivity
  Double_t EvalCountSensHalfLife(TString inFileName, TString varName, Double_t yrs, Double_t mass);
  // Discovery
  Double_t GetCriticalCounts(Double_t alpha, Double_t nBkgd);
  Double_t GetLeastDetectableSignal(Double_t prob, Double_t nBkgd, Double_t nCrit);
  Double_t EvalCountDiscHalfLife(TString inFileName, TString varName, Double_t yrs, Double_t mass);

  Double_t* ReadSensFromFiles(size_t n, double* yrs, const char* files, double mass);
  Double_t* ReadDiscFromFiles(size_t ny, double* yrs, size_t nc, double* cts, const char* files, double mass, Double_t prob, Double_t sigma, double max = 450., double min = -50., int bins = 10000, const char* treename = "tree", const char* cut = "stat_sig==0 && stat_bkg==0 && covQual_sig==3 && covQual_bkg==3");
  Double_t* EvalCountSensFromFiles(size_t n, double* yrs, const char* files, TString varName, double mass, double eff = 1.);
  Double_t* EvalCountDiscFromFiles(size_t n, double* yrs, const char* files, TString varName, double mass, double eff = 1.);
  
};

namespace nEXOPlotUtils
{
  
  // TCanvas* CreateHalflifeVsTimeCanvas(nEXONuclearMatrixElement::NME_t matElem, const char* name = "cc_[NME]", Double_t xMin=0., Double_t xMax=10., Double_t yMin=1.e25, Double_t yMax=2.e29);

}

#endif
