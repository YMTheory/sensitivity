#include "nEXOFitResult.hh"

ClassImp(nEXOFitResult)

nEXOFitResult::nEXOFitResult(const char* signalName) : fitres_sig(0), fitres_bkg(0)
{
  signal_name = signalName; // name of the signal variable (should be bb0n or LXeBb0n or something like that...)
 
}

nEXOFitResult::~nEXOFitResult()
{
  DeleteAllPointers();
}

void nEXOFitResult::DeleteAllPointers()
{
  if(fitres_sig)
    delete fitres_sig;
  fitres_sig = 0; // signal model fit result
  if(fitres_bkg)
    delete fitres_bkg;
  fitres_bkg = 0; //bkgd model fit result  
}

void nEXOFitResult::Reset()
{
  DeleteAllPointers();
  
  num_signal = 0; // best fit number of signal events
  num_signal_eHi = 0; // best fit upper error according to MINOS
  num_signal_eLo = 0; // best fit lower error according to MINOS 
  nll_sig = 0; // signal model NLL
  nll_bkg = 0; // bkgd model NLL
  nll_ratio = 0; // NLL ratio (signal/bkdg model)
  nll_offset = 0; //NLL offset
  stat_sig = 0; // migrad status for signal model fit
  covQual_sig = 0; // covariance quality code for signal model fit
  stat_bkg = 0; // migrad status for bkgd model fit
  covQual_bkg = 0; // covariance quality code for bkgd model fit
  real_time = 0; //real time to complete one fitting loop
  
  all_bkg = 0; // counts in SS bkg
  roi_bkg = 0; // counts in FWHM bkg
  roi_bkg_3t = 0; // counts in FWHM-3t bkg
  roi_bkg_1t = 0; // counts in FWHM-1t bkg
  
  bkg_tot = 0; // total bkgd before rebinning of PDFs
  bkg_fwhm_fv = 0; // FWHM bkg before rebinning of PDFs
  bkg_fwhm_3t = 0; // FWHM-1t bkg before rebinning of PDFs
  bkg_fwhm_1t = 0; // FWHM-3t bkg before rebinning of PDFs
}
