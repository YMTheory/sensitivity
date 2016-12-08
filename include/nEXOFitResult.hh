#ifndef nEXOFitResult_hh
#define nEXOFitResult_hh

#include <iostream>
#include <vector>

#include "TObject.h"
#include "TString.h"
#include "RooFitResult.h"


class nEXOFitResult : public TObject
{
public:

  nEXOFitResult(const char* signalName = 0);
  virtual ~nEXOFitResult();

  void Reset();

  TString signal_name; // name of the signal variable (should be bb0n or LXeBb0n or something like that...)
  Double_t num_signal; // best fit number of signal events
  Double_t num_signal_eHi; // best fit upper error according to MINOS
  Double_t num_signal_eLo; // best fit lower error according to MINOS 
  Double_t nll_sig; // signal model NLL
  Double_t nll_bkg; // bkgd model NLL
  Double_t nll_ratio; // NLL ratio (signal/bkdg model)
  Double_t nll_offset; //NLL offset
  Int_t stat_sig; // migrad status for signal model fit
  Int_t covQual_sig; // covariance quality code for signal model fit
  Int_t stat_bkg; // migrad status for bkgd model fit
  Int_t covQual_bkg; // covariance quality code for bkgd model fit
  RooFitResult *fitres_sig; // signal model fit result
  RooFitResult *fitres_bkg; //bkgd model fit result
  Double_t real_time; //real time to complete one fitting loop
  
  Double_t all_bkg; // counts in SS bkg after rebinning (not precisely the volume/energy region specified)
  Double_t roi_bkg; // counts in FWHM bkg after rebinning (not precisely the volume/energy region specified)
  Double_t roi_bkg_3t; // counts in FWHM-3t bkg after rebinning (not precisely the volume/energy region specified)
  Double_t roi_bkg_1t; // counts in FWHM-1t bkg after rebinning (not precisely the volume/energy region specified)

  Double_t bkg_tot; // total of mean bkgd per year before rebinning of PDFs
  Double_t bkg_fwhm_fv; // mean bkg/yr in FWHM before rebinning of PDFs
  Double_t bkg_fwhm_3t; // mean bkg/yr in FWHM-3t before rebinning of PDFs
  Double_t bkg_fwhm_2t; // mean bkg/yr in FWHM-2t before rebinning of PDFs
  Double_t bkg_fwhm_1t; // mean bkg/yr in FWHM-1t before rebinning of PDFs
  Double_t bkg_fwhm_3p5t; // mean bkg/yr in FWHM-3p5t before rebinning of PDFs
  Double_t bkg_fwhm_2p5t; // mean bkg/yr in FWHM-2p5t before rebinning of PDFs
  Double_t bkg_fwhm_1p5t; // mean bkg/yr in FWHM-1p5t before rebinning of PDFs
  Double_t bkg_fwhm_0p5t; // mean bkg/yr in FWHM-0p5t before rebinning of PDFs

  Double_t bkg_1sigma_fv; // mean bkg/yr in FWHM before rebinning of PDFs
  Double_t bkg_1sigma_3t; // mean bkg/yr in FWHM-3t before rebinning of PDFs
  Double_t bkg_1sigma_2t; // mean bkg/yr in FWHM-2t before rebinning of PDFs
  Double_t bkg_1sigma_1t; // mean bkg/yr in FWHM-1t before rebinning of PDFs
  Double_t bkg_1sigma_3p5t; // mean bkg/yr in FWHM-3p5t before rebinning of PDFs
  Double_t bkg_1sigma_2p5t; // mean bkg/yr in FWHM-2p5t before rebinning of PDFs
  Double_t bkg_1sigma_1p5t; // mean bkg/yr in FWHM-1p5t before rebinning of PDFs
  Double_t bkg_1sigma_0p5t; // mean bkg/yr in FWHM-0p5t before rebinning of PDFs

  Double_t bkg_2sigma_fv; // mean bkg/yr in FWHM before rebinning of PDFs
  Double_t bkg_2sigma_3t; // mean bkg/yr in FWHM-3t before rebinning of PDFs
  Double_t bkg_2sigma_2t; // mean bkg/yr in FWHM-2t before rebinning of PDFs
  Double_t bkg_2sigma_1t; // mean bkg/yr in FWHM-1t before rebinning of PDFs
  Double_t bkg_2sigma_3p5t; // mean bkg/yr in FWHM-3p5t before rebinning of PDFs
  Double_t bkg_2sigma_2p5t; // mean bkg/yr in FWHM-2p5t before rebinning of PDFs
  Double_t bkg_2sigma_1p5t; // mean bkg/yr in FWHM-1p5t before rebinning of PDFs
  Double_t bkg_2sigma_0p5t; // mean bkg/yr in FWHM-0p5t before rebinning of PDFs

  Double_t bkg_3sigma_fv; // mean bkg/yr in FWHM before rebinning of PDFs
  Double_t bkg_3sigma_3t; // mean bkg/yr in FWHM-3t before rebinning of PDFs
  Double_t bkg_3sigma_2t; // mean bkg/yr in FWHM-2t before rebinning of PDFs
  Double_t bkg_3sigma_1t; // mean bkg/yr in FWHM-1t before rebinning of PDFs
  Double_t bkg_3sigma_3p5t; // mean bkg/yr in FWHM-3p5t before rebinning of PDFs
  Double_t bkg_3sigma_2p5t; // mean bkg/yr in FWHM-2p5t before rebinning of PDFs
  Double_t bkg_3sigma_1p5t; // mean bkg/yr in FWHM-1p5t before rebinning of PDFs
  Double_t bkg_3sigma_0p5t; // mean bkg/yr in FWHM-0p5t before rebinning of PDFs

protected:
  void DeleteAllPointers();    
  
  ClassDef(nEXOFitResult,1)
};

#endif


