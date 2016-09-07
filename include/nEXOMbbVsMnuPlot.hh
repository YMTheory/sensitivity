#ifndef nEXOMbbVsMnuPlot_hh
#define nEXOMbbVsMnuPlot_hh

#include <iostream>
#include <map>

#include "TNamed.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGaxis.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"

#include "nEXOSensPlot.hh"
#include "nEXONuclearMatrixElement.hh"
#include "nEXOUtils.hh"

class nEXOMbbVsMnuPlot : public nEXOSensPlot
{
public:
  enum NuMassOption_t {kNuMassMin=0, kNuMassSum};
  static const size_t nNuMassOptions = 2;

  nEXOMbbVsMnuPlot(const char* name = 0, const char* title = 0);
  virtual ~nEXOMbbVsMnuPlot();

  bool SetBand(const char* name, const char* title, double min, double max);  
  bool SetBandFromFile(const char* name, const char* title, const char* filename, double yrs, double mass);

  TObject* GetPlot();

  TCanvas* CreateEmptyCanvas();
  
  void PlotEXO200(TCanvas& canvas);
  
  std::map<TString, Int_t> fLineColors;
  std::map<TString, Int_t> fLineWidths;
  std::map<TString, Int_t> fFillColors;
  std::map<TString, Int_t> fFillStyles;
  std::map<TString, Int_t> fTextColors;
  std::map<TString, Int_t> fTextFonts;
  std::map<TString, Double_t> fTextSizes;

  NuMassOption_t fNuMassOption;

protected:
  
  Int_t fNpoints;
  Double_t fLoXRangeCoeff;
  Int_t fLoXRangeExp; // assumed negative (power)
  Double_t fHiXRangeCoeff;
  Int_t fHiXRangeExp; // assumed positive (power)
  Double_t fLoYRangeCoeff;
  Int_t fLoYRangeExp; // assumed negative (power)
  Double_t fHiYRangeCoeff;
  Int_t fHiYRangeExp; // assumed positive (power)

  std::map<TString, TGraph> fBands;
  std::map<TString, TPaveText> fTitles;

  Double_t GetLoXRange(){return fLoXRangeCoeff*TMath::Power(10., -fLoXRangeExp);}
  Double_t GetHiXRange(){return fHiXRangeCoeff*TMath::Power(10., fHiXRangeExp);};
  Double_t GetLoYRange(){return fLoYRangeCoeff*TMath::Power(10., -fLoYRangeExp);};
  Double_t GetHiYRange(){return fHiYRangeCoeff*TMath::Power(10., fHiYRangeExp);};
  void CreateAvailableBands();
  void FillGraphPointsMmin(TGraph* g_n_cv,TGraph* g_n_err,TGraph* g_i_cv,TGraph* g_i_err);
  void FillGraphPointsMsum(TGraph* g_n_cv,TGraph* g_n_err,TGraph* g_i_cv,TGraph* g_i_err);

  ClassDef(nEXOMbbVsMnuPlot,1) // Class for mbb vs mnu plots
};

#endif

