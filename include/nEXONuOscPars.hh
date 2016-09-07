#ifndef nEXONuOscPars_hh
#define nEXONuOscPars_hh

#include <iostream>

#include "TROOT.h"
#include "TF3.h"
#include "Math/MultiRootFinder.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/WrappedTF1.h"
#include "Math/RootFinderAlgorithms.h"
#include "Math/RootFinder.h"

class nEXONuOscPars
{
public:
  enum NuOscParsFit_t { k1205_4018v4=0 , k1205_5254 , kprd90_093006_2014};

  static nEXONuOscPars* GetInstance(){
    if(nEXONuOscPars::fInstance == NULL) nEXONuOscPars::fInstance = new nEXONuOscPars();
    return nEXONuOscPars::fInstance;
  }

  virtual ~nEXONuOscPars(){};

  bool SetFitSource(NuOscParsFit_t fitSource);
  void EvalNormalMbbMmin(Double_t* mbb, Double_t* mbb_lo, Double_t* mbb_hi, Double_t m_min_v, Int_t nCL);
  void EvalInvertedMbbMmin(Double_t* mbb, Double_t* mbb_lo, Double_t* mbb_hi, Double_t m_min_v, Int_t nCL);
  bool EvalNormalMbbMsum(Double_t* mbb, Double_t* mbb_lo, Double_t* mbb_hi, Double_t m_sum_v, Int_t nCL);
  bool EvalInvertedMbbMsum(Double_t* mbb, Double_t* mbb_lo, Double_t* mbb_hi, Double_t m_sum_v, Int_t nCL);

  void SetNormalNuPars(Double_t& s12,Double_t& s13,Double_t& m12,Double_t& m13,Double_t& s12_lo,Double_t& s13_lo,Double_t& m12_lo,Double_t& m13_lo,Double_t& s12_hi,Double_t& s13_hi,Double_t& m12_hi,Double_t& m13_hi, Int_t nCL);
  void SetInvertedNuPars(Double_t& s12,Double_t& s13,Double_t& m12,Double_t& m13,Double_t& s12_lo,Double_t& s13_lo,Double_t& m12_lo,Double_t& m13_lo,Double_t& s12_hi,Double_t& s13_hi,Double_t& m12_hi,Double_t& m13_hi, Int_t nCL);
  void EvalNormalMassesMmin(Double_t& m1,Double_t& m2,Double_t& m3,Double_t m12,Double_t m13,Double_t m_min_v);
  void EvalInvertedMassesMmin(Double_t& m1,Double_t& m2,Double_t& m3,Double_t m12,Double_t m13,Double_t m_min_v);
  bool EvalNormalMassesMsum(Double_t& m1,Double_t& m2,Double_t& m3,Double_t m12,Double_t m13,Double_t m_sum_v);
  bool EvalInvertedMassesMsum(Double_t& m1,Double_t& m2,Double_t& m3,Double_t m12,Double_t m13,Double_t m_sum_v);

private:
  nEXONuOscPars(NuOscParsFit_t fitSource = k1205_4018v4);
  static nEXONuOscPars* fInstance;

  NuOscParsFit_t fFitSource;

  Double_t s12_cv;	Double_t m12_cv;
  Double_t s13_n_cv;	Double_t m13_n_cv;
  Double_t s13_i_cv;	Double_t m13_i_cv;
  
  Double_t s12_68[2];	Double_t m12_68[2];
  Double_t s13_n_68[2];	Double_t m13_n_68[2];
  Double_t s13_i_68[2];	Double_t m13_i_68[2];
  
  Double_t s12_90[2];	Double_t m12_90[2];
  Double_t s13_n_90[2];	Double_t m13_n_90[2];
  Double_t s13_i_90[2];	Double_t m13_i_90[2];
  
  Double_t s12_95[2];	Double_t m12_95[2];
  Double_t s13_n_95[2];	Double_t m13_n_95[2];
  Double_t s13_i_95[2];	Double_t m13_i_95[2];

  bool SetFitSource_1205_4018v4();
  bool SetFitSource_1205_5254();
  bool SetFitSource_prd90_093006_2014();

  ClassDef(nEXONuOscPars,1)	// Class to carry neutrino oscillation parameters
};

#endif
