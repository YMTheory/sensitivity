#include "nEXONuOscPars.hh"

ClassImp(nEXONuOscPars);

nEXONuOscPars* nEXONuOscPars::fInstance = 0;

nEXONuOscPars::nEXONuOscPars(NuOscParsFit_t fitSource)
{
  SetFitSource(fitSource);
}

bool nEXONuOscPars::SetFitSource(NuOscParsFit_t fitSource)
{
  switch(fitSource)
  {
      case k1205_4018v4:
        return SetFitSource_1205_4018v4();
      case k1205_5254:
        return SetFitSource_1205_5254();
      case kprd90_093006_2014:
        return SetFitSource_prd90_093006_2014();
  }
  return false;
}

bool nEXONuOscPars::SetFitSource_1205_4018v4()
{
  //Global fits to neutrino oscillation parameters from arXiv:1205.4018v4
  
  s12_cv 	= 0.320;		m12_cv 	= 7.62e-5;
  s13_n_cv = 0.0246;	m13_n_cv	= 2.55e-3;
  s13_i_cv = 0.0250;	m13_i_cv	= 2.43e-3;
  
  s12_68	[0] = 0.3030; 	s12_68	[1] = 0.3360;
  m12_68	[0] = 7.43e-5; m12_68	[1] = 7.81e-5;
  s13_n_68	[0] = 0.0218; 	s13_n_68	[1] = 0.0275;
  m13_n_68	[0] = 2.46e-3; m13_n_68	[1] = 2.61e-3;
  s13_i_68	[0] = 0.0223;	s13_i_68	[1] = 0.0276;
  m13_i_68	[0] = 2.37e-3; m13_i_68	[1] = 2.50e-3;
  
  s12_90	[0] = 0.2900;	s12_90	[1] = 0.3500;
  m12_90	[0] = 7.27e-5; m12_90	[1] = 8.01e-5;
  s13_n_90	[0] = 0.0190;	s13_n_90	[1] = 0.0300;
  m13_n_90	[0] = 2.38e-3; m13_n_90	[1] = 2.68e-3;
  s13_i_90	[0] = 0.0200;	s13_i_90	[1] = 0.0300;
  m13_i_90	[0] = 2.29e-3; m13_i_90	[1] = 2.58e-3;
  
  s12_95	[0] = 0.2700;	s12_95	[1] = 0.3700;
  m12_95	[0] = 7.12e-5;	m12_95	[1] = 8.20e-5;
  s13_n_95	[0] = 0.0170;	s13_n_95	[1] = 0.0330;
  m13_n_95	[0] = 2.31e-3;	m13_n_95	[1] = 2.74e-3;
  s13_i_95	[0] = 0.0170; 	s13_i_95	[1] = 0.0330;
  m13_i_95	[0] = 2.21e-3; m13_i_95	[1] = 2.64e-3;

  return true;
}

bool nEXONuOscPars::SetFitSource_1205_5254()
{
  //Global fits to neutrino oscillation parameters from arXiv:1205.5254

  s12_cv 	= 0.307;		m12_cv 	= 7.54e-5;
  s13_n_cv = 0.0241;	m13_n_cv	= 2.43e-3;
  s13_i_cv = 0.0244;	m13_i_cv	= 2.42e-3;
  
  s12_68	[0] = 0.291; 	s12_68	[1] = 0.325;
  m12_68	[0] = 7.32e-5; m12_68	[1] = 7.80e-5;
  s13_n_68	[0] = 0.0216; 	s13_n_68	[1] = 0.0266;
  m13_n_68	[0] = 2.33e-3; m13_n_68	[1] = 2.49e-3;
  s13_i_68	[0] = 0.0219;	s13_i_68	[1] = 0.0267;
  m13_i_68	[0] = 2.31e-3; m13_i_68	[1] = 2.49e-3;
  
  s12_90	[0] = 0.275;	s12_90	[1] = 0.3420;
  m12_90	[0] = 7.15e-5; m12_90	[1] = 8.00e-5;
  s13_n_90	[0] = 0.0193;	s13_n_90	[1] = 0.0290;
  m13_n_90	[0] = 2.38e-3; m13_n_90	[1] = 2.68e-3;
  s13_i_90	[0] = 0.0194;	s13_i_90	[1] = 0.0291;
  m13_i_90	[0] = 2.26e-3; m13_i_90	[1] = 2.53e-3;
  
  s12_95	[0] = 0.2700;	s12_95	[1] = 0.3700;
  m12_95	[0] = 7.12e-5;	m12_95	[1] = 8.20e-5;
  s13_n_95	[0] = 0.0170;	s13_n_95	[1] = 0.0330;
  m13_n_95	[0] = 2.31e-3;	m13_n_95	[1] = 2.74e-3;
  s13_i_95	[0] = 0.0170; 	s13_i_95	[1] = 0.0330;
  m13_i_95	[0] = 2.21e-3; m13_i_95	[1] = 2.64e-3;

  return true;
}

bool nEXONuOscPars::SetFitSource_prd90_093006_2014()
{
  //Global fits to neutrino oscillation parameters from PRD 90 093006 (2014)
  s12_cv = 0.323;		m12_cv 	= 7.60e-5;
  s13_n_cv = 0.0234;	m13_n_cv	= 2.48e-3;
  s13_i_cv = 0.0240;	m13_i_cv	= 2.38e-3;
  
  s12_68	[0] = 0.307; 	s12_68	[1] = 0.339;
  m12_68	[0] = 7.42e-5; m12_68	[1] = 7.79e-5;
  s13_n_68	[0] = 0.0214; 	s13_n_68[1] = 0.0254;
  m13_n_68	[0] = 2.41e-3; m13_n_68	[1] = 2.53e-3;
  s13_i_68	[0] = 0.0221;	s13_i_68[1] = 0.0259;
  m13_i_68	[0] = 2.32e-3; m13_i_68	[1] = 2.43e-3;
  
  s12_90	[0] = 0.292;	s12_90	[1] = 0.357;
  m12_90	[0] = 7.26e-5; m12_90	[1] = 7.99e-5;
  s13_n_90	[0] = 0.0195;	s13_n_90[1] = 0.0274;
  m13_n_90	[0] = 2.35e-3; m13_n_90	[1] = 2.59e-3;
  s13_i_90	[0] = 0.0202;	s13_i_90[1] = 0.0278;
  m13_i_90	[0] = 2.26e-3; m13_i_90	[1] = 2.48e-3;
  
  s12_95	[0] = 0.278;	s12_95	[1] = 0.375;
  m12_95	[0] = 7.11e-5;	m12_95	[1] = 8.18e-5;
  s13_n_95	[0] = 0.0177;	s13_n_95[1] = 0.0294;
  m13_n_95	[0] = 2.30e-3;	m13_n_95[1] = 2.65e-3;
  s13_i_95	[0] = 0.0183; 	s13_i_95[1] = 0.0297;
  m13_i_95	[0] = 2.20e-3; m13_i_95	[1] = 2.54e-3;

  return true;
}

void nEXONuOscPars::EvalNormalMbbMmin(Double_t* mbb, Double_t* mbb_lo, Double_t* mbb_hi, Double_t m_min_v, Int_t nCL)
{
  
  Double_t s12=0., s13=0., m12=0., m13=0.;
  Double_t s12_lo=0., s13_lo=0., m12_lo=0., m13_lo=0.;
  Double_t s12_hi=0., s13_hi=0., m12_hi=0., m13_hi=0.;

  SetNormalNuPars(s12,s13,m12,m13,s12_lo,s13_lo,m12_lo,m13_lo,s12_hi,s13_hi,m12_hi,m13_hi,nCL);

  Double_t m1 = 0.;
  Double_t m2 = 0.;
  Double_t m3 = 0.;
  
  EvalNormalMassesMmin(m1,m2,m3,m12,m13,m_min_v);

  Double_t c12 = 1. - s12;
  Double_t c13 = 1. - s13;
  
  mbb[0] = TMath::Abs(c12*c13*m1 + s12*c13*m2 + s13*m3);
  mbb[1] = TMath::Abs(c12*c13*m1 + s12*c13*m2 - s13*m3);
  mbb[2] = TMath::Abs(c12*c13*m1 - s12*c13*m2 + s13*m3);
  mbb[3] = TMath::Abs(c12*c13*m1 - s12*c13*m2 - s13*m3);
  
  Double_t term_s12_lo_p = c13*c13*(m1-m2)*(m1-m2)*s12_lo*s12_lo;
  Double_t term_s13_lo_p = (c12*m1+s12*m2)*(c12*m1+s12*m2)*s13_lo*s13_lo;
  Double_t term_s12_lo_m = c13*c13*(m1+m2)*(m1+m2)*s12_lo*s12_lo;
  Double_t term_s13_lo_m = (c12*m1-s12*m2)*(c12*m1-s12*m2)*s13_lo*s13_lo;
  Double_t term_m12_lo = s12*s12*c13*c13/(4.*m2*m2)*m12_lo*m12_lo;
  Double_t term_m13_lo = s13*s13/(4.*m3*m3)*m13_lo*m13_lo;
  mbb_lo[0] = sqrt(term_s12_lo_p + term_s13_lo_p + term_m12_lo + term_m13_lo);
  mbb_lo[1] = sqrt(term_s12_lo_p + term_s13_lo_p + term_m12_lo + term_m13_lo);
  mbb_lo[2] = sqrt(term_s12_lo_m + term_s13_lo_m + term_m12_lo + term_m13_lo);
  mbb_lo[3] = sqrt(term_s12_lo_m + term_s13_lo_m + term_m12_lo + term_m13_lo);
  
  Double_t term_s12_hi_p = c13*c13*(m1-m2)*(m1-m2)*s12_hi*s12_hi;
  Double_t term_s13_hi_p = (c12*m1+s12*m2)*(c12*m1+s12*m2)*s13_hi*s13_hi;
  Double_t term_s12_hi_m = c13*c13*(m1+m2)*(m1+m2)*s12_hi*s12_hi;
  Double_t term_s13_hi_m = (c12*m1-s12*m2)*(c12*m1-s12*m2)*s13_hi*s13_hi;
  Double_t term_m12_hi = s12*s12*c13*c13/(4.*m2*m2)*m12_hi*m12_hi;
  Double_t term_m13_hi = s13*s13/(4.*m3*m3)*m13_hi*m13_hi;
  mbb_hi[0] = sqrt(term_s12_hi_p + term_s13_hi_p + term_m12_hi + term_m13_hi);
  mbb_hi[1] = sqrt(term_s12_hi_p + term_s13_hi_p + term_m12_hi + term_m13_hi);
  mbb_hi[2] = sqrt(term_s12_hi_m + term_s13_hi_m + term_m12_hi + term_m13_hi);
  mbb_hi[3] = sqrt(term_s12_hi_m + term_s13_hi_m + term_m12_hi + term_m13_hi);
}

void nEXONuOscPars::EvalInvertedMbbMmin(Double_t* mbb, Double_t* mbb_lo, Double_t* mbb_hi, Double_t m_min_v, Int_t nCL)
{
  Double_t s12=0., s13=0., m12=0., m13=0.;
  Double_t s12_lo=0., s13_lo=0., m12_lo=0., m13_lo=0.;
  Double_t s12_hi=0., s13_hi=0., m12_hi=0., m13_hi=0.;

  SetInvertedNuPars(s12,s13,m12,m13,s12_lo,s13_lo,m12_lo,m13_lo,s12_hi,s13_hi,m12_hi,m13_hi,nCL);

  Double_t m1 = 0.;
  Double_t m2 = 0.;
  Double_t m3 = 0.;
  
  EvalInvertedMassesMmin(m1,m2,m3,m12,m13,m_min_v);
    
  Double_t c12 = 1. - s12;
  Double_t c13 = 1. - s13;  

  mbb[0] = TMath::Abs(c12*c13*m1 + s12*c13*m2 + s13*m3);
  mbb[1] = TMath::Abs(c12*c13*m1 + s12*c13*m2 - s13*m3);
  mbb[2] = TMath::Abs(c12*c13*m1 - s12*c13*m2 + s13*m3);
  mbb[3] = TMath::Abs(c12*c13*m1 - s12*c13*m2 - s13*m3);
  
  Double_t term_s12_lo_p = c13*c13*(m1-m2)*(m1-m2)*s12_lo*s12_lo;
  Double_t term_s13_lo_p = (c12*m1+s12*m2)*(c12*m1+s12*m2)*s13_lo*s13_lo;
  Double_t term_s12_lo_m = c13*c13*(m1+m2)*(m1+m2)*s12_lo*s12_lo;
  Double_t term_s13_lo_m = (c12*m1-s12*m2)*(c12*m1-s12*m2)*s13_lo*s13_lo;
  Double_t term_m12_lo = s12*s12*c13*c13/(4.*m2*m2)*m12_lo*m12_lo;
  Double_t term_m13_lo_p = (s12*c13/2./m1 + s12*c13/2./m2)*(s12*c13/2./m1 + s12*c13/2./m2)*m13_lo*m13_lo;
  Double_t term_m13_lo_m = (s12*c13/2./m1 - s12*c13/2./m2)*(s12*c13/2./m1 - s12*c13/2./m2)*m13_lo*m13_lo;
  mbb_lo[0] = sqrt(term_s12_lo_p + term_s13_lo_p + term_m12_lo + term_m13_lo_p);
  mbb_lo[1] = sqrt(term_s12_lo_p + term_s13_lo_p + term_m12_lo + term_m13_lo_p);
  mbb_lo[2] = sqrt(term_s12_lo_m + term_s13_lo_m + term_m12_lo + term_m13_lo_m);
  mbb_lo[3] = sqrt(term_s12_lo_m + term_s13_lo_m + term_m12_lo + term_m13_lo_m);
  
  Double_t term_s12_hi_p = c13*c13*(m1-m2)*(m1-m2)*s12_hi*s12_hi;
  Double_t term_s13_hi_p = (c12*m1+s12*m2)*(c12*m1+s12*m2)*s13_hi*s13_hi;
  Double_t term_s12_hi_m = c13*c13*(m1+m2)*(m1+m2)*s12_hi*s12_hi;
  Double_t term_s13_hi_m = (c12*m1-s12*m2)*(c12*m1-s12*m2)*s13_hi*s13_hi;
  Double_t term_m12_hi = s12*s12*c13*c13/(4.*m2*m2)*m12_hi*m12_hi;
  Double_t term_m13_hi_p = (s12*c13/2./m1 + s12*c13/2./m2)*(s12*c13/2./m1 + s12*c13/2./m2)*m13_hi*m13_hi;
  Double_t term_m13_hi_m = (s12*c13/2./m1 - s12*c13/2./m2)*(s12*c13/2./m1 - s12*c13/2./m2)*m13_hi*m13_hi;
  mbb_hi[0] = sqrt(term_s12_hi_p + term_s13_hi_p + term_m12_hi + term_m13_hi_p);
  mbb_hi[1] = sqrt(term_s12_hi_p + term_s13_hi_p + term_m12_hi + term_m13_hi_p);
  mbb_hi[2] = sqrt(term_s12_hi_m + term_s13_hi_m + term_m12_hi + term_m13_hi_m);
  mbb_hi[3] = sqrt(term_s12_hi_m + term_s13_hi_m + term_m12_hi + term_m13_hi_m);
}

bool nEXONuOscPars::EvalNormalMbbMsum(Double_t* mbb, Double_t* mbb_lo, Double_t* mbb_hi, Double_t m_sum_v, Int_t nCL)
{
  
  Double_t s12=0., s13=0., m12=0., m13=0.;
  Double_t s12_lo=0., s13_lo=0., m12_lo=0., m13_lo=0.;
  Double_t s12_hi=0., s13_hi=0., m12_hi=0., m13_hi=0.;

  SetNormalNuPars(s12,s13,m12,m13,s12_lo,s13_lo,m12_lo,m13_lo,s12_hi,s13_hi,m12_hi,m13_hi,nCL);

  Double_t m1 = 0.;
  Double_t m2 = 0.;
  Double_t m3 = 0.;
  
  if(not EvalNormalMassesMsum(m1,m2,m3,m12,m13,m_sum_v))
    return false;

  Double_t c12 = 1. - s12;
  Double_t c13 = 1. - s13;
  
  mbb[0] = TMath::Abs(c12*c13*m1 + s12*c13*m2 + s13*m3);
  mbb[1] = TMath::Abs(c12*c13*m1 + s12*c13*m2 - s13*m3);
  mbb[2] = TMath::Abs(c12*c13*m1 - s12*c13*m2 + s13*m3);
  mbb[3] = TMath::Abs(c12*c13*m1 - s12*c13*m2 - s13*m3);
  
  Double_t term_s12_lo_p = c13*c13*(m1-m2)*(m1-m2)*s12_lo*s12_lo;
  Double_t term_s13_lo_p = (c12*m1+s12*m2)*(c12*m1+s12*m2)*s13_lo*s13_lo;
  Double_t term_s12_lo_m = c13*c13*(m1+m2)*(m1+m2)*s12_lo*s12_lo;
  Double_t term_s13_lo_m = (c12*m1-s12*m2)*(c12*m1-s12*m2)*s13_lo*s13_lo;
  Double_t term_m12_lo = s12*s12*c13*c13/(4.*m2*m2)*m12_lo*m12_lo;
  Double_t term_m13_lo = s13*s13/(4.*m3*m3)*m13_lo*m13_lo;
  mbb_lo[0] = sqrt(term_s12_lo_p + term_s13_lo_p + term_m12_lo + term_m13_lo);
  mbb_lo[1] = sqrt(term_s12_lo_p + term_s13_lo_p + term_m12_lo + term_m13_lo);
  mbb_lo[2] = sqrt(term_s12_lo_m + term_s13_lo_m + term_m12_lo + term_m13_lo);
  mbb_lo[3] = sqrt(term_s12_lo_m + term_s13_lo_m + term_m12_lo + term_m13_lo);
  
  Double_t term_s12_hi_p = c13*c13*(m1-m2)*(m1-m2)*s12_hi*s12_hi;
  Double_t term_s13_hi_p = (c12*m1+s12*m2)*(c12*m1+s12*m2)*s13_hi*s13_hi;
  Double_t term_s12_hi_m = c13*c13*(m1+m2)*(m1+m2)*s12_hi*s12_hi;
  Double_t term_s13_hi_m = (c12*m1-s12*m2)*(c12*m1-s12*m2)*s13_hi*s13_hi;
  Double_t term_m12_hi = s12*s12*c13*c13/(4.*m2*m2)*m12_hi*m12_hi;
  Double_t term_m13_hi = s13*s13/(4.*m3*m3)*m13_hi*m13_hi;
  mbb_hi[0] = sqrt(term_s12_hi_p + term_s13_hi_p + term_m12_hi + term_m13_hi);
  mbb_hi[1] = sqrt(term_s12_hi_p + term_s13_hi_p + term_m12_hi + term_m13_hi);
  mbb_hi[2] = sqrt(term_s12_hi_m + term_s13_hi_m + term_m12_hi + term_m13_hi);
  mbb_hi[3] = sqrt(term_s12_hi_m + term_s13_hi_m + term_m12_hi + term_m13_hi);

  return true;
}

bool nEXONuOscPars::EvalInvertedMbbMsum(Double_t* mbb, Double_t* mbb_lo, Double_t* mbb_hi, Double_t m_sum_v, Int_t nCL)
{
  Double_t s12=0., s13=0., m12=0., m13=0.;
  Double_t s12_lo=0., s13_lo=0., m12_lo=0., m13_lo=0.;
  Double_t s12_hi=0., s13_hi=0., m12_hi=0., m13_hi=0.;

  SetInvertedNuPars(s12,s13,m12,m13,s12_lo,s13_lo,m12_lo,m13_lo,s12_hi,s13_hi,m12_hi,m13_hi,nCL);

  Double_t m1 = 0.;
  Double_t m2 = 0.;
  Double_t m3 = 0.;
  
  if(not EvalInvertedMassesMsum(m1,m2,m3,m12,m13,m_sum_v))
    return false;
    
  Double_t c12 = 1. - s12;
  Double_t c13 = 1. - s13;  

  mbb[0] = TMath::Abs(c12*c13*m1 + s12*c13*m2 + s13*m3);
  mbb[1] = TMath::Abs(c12*c13*m1 + s12*c13*m2 - s13*m3);
  mbb[2] = TMath::Abs(c12*c13*m1 - s12*c13*m2 + s13*m3);
  mbb[3] = TMath::Abs(c12*c13*m1 - s12*c13*m2 - s13*m3);
  
  Double_t term_s12_lo_p = c13*c13*(m1-m2)*(m1-m2)*s12_lo*s12_lo;
  Double_t term_s13_lo_p = (c12*m1+s12*m2)*(c12*m1+s12*m2)*s13_lo*s13_lo;
  Double_t term_s12_lo_m = c13*c13*(m1+m2)*(m1+m2)*s12_lo*s12_lo;
  Double_t term_s13_lo_m = (c12*m1-s12*m2)*(c12*m1-s12*m2)*s13_lo*s13_lo;
  Double_t term_m12_lo = s12*s12*c13*c13/(4.*m2*m2)*m12_lo*m12_lo;
  Double_t term_m13_lo_p = (s12*c13/2./m1 + s12*c13/2./m2)*(s12*c13/2./m1 + s12*c13/2./m2)*m13_lo*m13_lo;
  Double_t term_m13_lo_m = (s12*c13/2./m1 - s12*c13/2./m2)*(s12*c13/2./m1 - s12*c13/2./m2)*m13_lo*m13_lo;
  mbb_lo[0] = sqrt(term_s12_lo_p + term_s13_lo_p + term_m12_lo + term_m13_lo_p);
  mbb_lo[1] = sqrt(term_s12_lo_p + term_s13_lo_p + term_m12_lo + term_m13_lo_p);
  mbb_lo[2] = sqrt(term_s12_lo_m + term_s13_lo_m + term_m12_lo + term_m13_lo_m);
  mbb_lo[3] = sqrt(term_s12_lo_m + term_s13_lo_m + term_m12_lo + term_m13_lo_m);
  
  Double_t term_s12_hi_p = c13*c13*(m1-m2)*(m1-m2)*s12_hi*s12_hi;
  Double_t term_s13_hi_p = (c12*m1+s12*m2)*(c12*m1+s12*m2)*s13_hi*s13_hi;
  Double_t term_s12_hi_m = c13*c13*(m1+m2)*(m1+m2)*s12_hi*s12_hi;
  Double_t term_s13_hi_m = (c12*m1-s12*m2)*(c12*m1-s12*m2)*s13_hi*s13_hi;
  Double_t term_m12_hi = s12*s12*c13*c13/(4.*m2*m2)*m12_hi*m12_hi;
  Double_t term_m13_hi_p = (s12*c13/2./m1 + s12*c13/2./m2)*(s12*c13/2./m1 + s12*c13/2./m2)*m13_hi*m13_hi;
  Double_t term_m13_hi_m = (s12*c13/2./m1 - s12*c13/2./m2)*(s12*c13/2./m1 - s12*c13/2./m2)*m13_hi*m13_hi;
  mbb_hi[0] = sqrt(term_s12_hi_p + term_s13_hi_p + term_m12_hi + term_m13_hi_p);
  mbb_hi[1] = sqrt(term_s12_hi_p + term_s13_hi_p + term_m12_hi + term_m13_hi_p);
  mbb_hi[2] = sqrt(term_s12_hi_m + term_s13_hi_m + term_m12_hi + term_m13_hi_m);
  mbb_hi[3] = sqrt(term_s12_hi_m + term_s13_hi_m + term_m12_hi + term_m13_hi_m);

  return true;
}

void nEXONuOscPars::SetNormalNuPars(Double_t& s12,Double_t& s13,Double_t& m12,Double_t& m13,Double_t& s12_lo,Double_t& s13_lo,Double_t& m12_lo,Double_t& m13_lo,Double_t& s12_hi,Double_t& s13_hi,Double_t& m12_hi,Double_t& m13_hi, Int_t nCL)
{
  // set nu-oscillation parameters for normal ordering with nCL limits
  s12 = s12_cv; s13 = s13_n_cv; m12 = m12_cv; m13 = m13_n_cv;
  
  switch (nCL) {
      case 1 :
        s12_lo = -s12_68[0]+s12; s13_lo = -s13_n_68[0]+s13; m12_lo = -m12_68[0]+m12; m13_lo = -m13_n_68[0]+m13;
        s12_hi = +s12_68[1]-s12; s13_hi = +s13_n_68[1]-s13; m12_hi = +m12_68[1]-m12; m13_hi = +m13_n_68[1]-m13;
        break;
      case 2 :
        s12_lo = -s12_90[0]+s12; s13_lo = -s13_n_90[0]+s13; m12_lo = -m12_90[0]+m12; m13_lo = -m13_n_90[0]+m13;
        s12_hi = +s12_90[1]-s12; s13_hi = +s13_n_90[1]-s13; m12_hi = +m12_90[1]-m12; m13_hi = +m13_n_90[1]-m13;
        break;
      case 3 :
        s12_lo = -s12_95[0]+s12; s13_lo = -s13_n_95[0]+s13; m12_lo = -m12_95[0]+m12; m13_lo = -m13_n_95[0]+m13;
        s12_hi = +s12_95[1]-s12; s13_hi = +s13_n_95[1]-s13; m12_hi = +m12_95[1]-m12; m13_hi = +m13_n_95[1]-m13;
        break;
  }	
}

void nEXONuOscPars::SetInvertedNuPars(Double_t& s12,Double_t& s13,Double_t& m12,Double_t& m13,Double_t& s12_lo,Double_t& s13_lo,Double_t& m12_lo,Double_t& m13_lo,Double_t& s12_hi,Double_t& s13_hi,Double_t& m12_hi,Double_t& m13_hi, Int_t nCL)
{
  // set nu-oscillation parameters for inverted ordering with nCL limits
  s12 = s12_cv; s13 = s13_i_cv; m12 = m12_cv; m13 = m13_i_cv;
  
  switch (nCL) {
      case 1 :
        s12_lo = -s12_68[0]+s12; s13_lo = -s13_i_68[0]+s13; m12_lo = -m12_68[0]+m12; m13_lo = -m13_i_68[0]+m13;
        s12_hi = +s12_68[1]-s12; s13_hi = +s13_i_68[1]-s13; m12_hi = +m12_68[1]-m12; m13_hi = +m13_i_68[1]-m13;
        break;
      case 2 :
        s12_lo = -s12_90[0]+s12; s13_lo = -s13_i_90[0]+s13; m12_lo = -m12_90[0]+m12; m13_lo = -m13_i_90[0]+m13;
        s12_hi = +s12_90[1]-s12; s13_hi = +s13_i_90[1]-s13; m12_hi = +m12_90[1]-m12; m13_hi = +m13_i_90[1]-m13;
        break;
      case 3 :
        s12_lo = -s12_95[0]+s12; s13_lo = -s13_i_95[0]+s13; m12_lo = -m12_95[0]+m12; m13_lo = -m13_i_95[0]+m13;
        s12_hi = +s12_95[1]-s12; s13_hi = +s13_i_95[1]-s13; m12_hi = +m12_95[1]-m12; m13_hi = +m13_i_95[1]-m13;
        break;
  }
  
}

void nEXONuOscPars::EvalNormalMassesMmin(Double_t& m1,Double_t& m2,Double_t& m3,Double_t m12,Double_t m13,Double_t m_min_v)
{
  m1 = m_min_v;
  m2 = sqrt(m1*m1 + m12);
  m3 = sqrt(m1*m1 + m13);
}

void nEXONuOscPars::EvalInvertedMassesMmin(Double_t& m1,Double_t& m2,Double_t& m3,Double_t m12,Double_t m13,Double_t m_min_v)
{
  m3 = m_min_v;
  m1 = sqrt(m3*m3 + m13);
  m2 = sqrt(m3*m3 + m12 + m13);
}

bool nEXONuOscPars::EvalNormalMassesMsum(Double_t& m1,Double_t& m2,Double_t& m3,Double_t m12,Double_t m13,Double_t m_sum_v)
{
  if(m_sum_v < sqrt(m12) + sqrt(m13))
    return false;

  TF1 fm1("fm1","[0] - x - sqrt([1] + x*x) - sqrt([2] + x*x)",0,1);
  fm1.SetParameters(m_sum_v,m12,m13);
  ROOT::Math::WrappedTF1 wfm1(fm1);

  ROOT::Math::RootFinder rf;//(ROOT::Math::RootFinder::kGSL_BISECTION);
  rf.SetFunction(wfm1,0,1);
  rf.Solve(10000);

  m1 = rf.Root();
  m2 = sqrt(m1*m1 + m12);
  m3 = sqrt(m1*m1 + m13);

  if(isnan(m1) or isnan(m2) or isnan(m3))
    return false;     

  return true;
}

bool nEXONuOscPars::EvalInvertedMassesMsum(Double_t& m1,Double_t& m2,Double_t& m3,Double_t m12,Double_t m13,Double_t m_sum_v)
{
  TF1 fm1("fm1","[0] - x - sqrt([1] + x*x) - sqrt(x*x - [2])",0,1);
  fm1.SetParameters(m_sum_v,m12,m13);
  ROOT::Math::WrappedTF1 wfm1(fm1);

  ROOT::Math::RootFinder rf;//(ROOT::Math::RootFinder::kGSL_BISECTION);
  rf.SetFunction(wfm1,0,1);
  rf.Solve(10000);

  m1 = rf.Root();
  m2 = sqrt(m1*m1 + m12);
  m3 = sqrt(m1*m1 - m13);

  if(isnan(m1) or isnan(m2) or isnan(m3))
    return false;     
  
  return true;
}
