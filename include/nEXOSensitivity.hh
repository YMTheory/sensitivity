#ifndef nEXOSensitivity_hh
#define nEXOSensitivity_hh

#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>

#include "TObject.h"
#include "TString.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TRint.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "Math/QuantFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooNLLVar.h"
#include "RooAddition.h"
#include "RooMinuit.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooMultiVarGaussian.h"
#include "RooFitResult.h"
#include "RooRandom.h"
#include "RooMsgService.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooConstVar.h"
#include "RooGenericPdf.h"
#include "RooBinning.h"

#include "nEXOFitResult.hh"

typedef struct ExcelTableValues {
  TString fPdf, fComponent, fIsotope, fMC_ID;
  Double_t fHalflife;

  std::vector<TString> fSuffixes;
  std::vector<Double_t> fCountsCV, fCountsError, fCountsUL;
  std::vector<Double_t> fRatioFWHMfv, fRatioFWHM3t, fRatioFWHM2t, fRatioFWHM1t, fRatioFWHM3p5t, fRatioFWHM2p5t, fRatioFWHM1p5t, fRatioFWHM0p5t;
  std::vector<Double_t> fRatio1sfv, fRatio1s3t, fRatio1s2t, fRatio1s1t, fRatio1s3p5t, fRatio1s2p5t, fRatio1s1p5t, fRatio1s0p5t;
  std::vector<Double_t> fRatio2sfv, fRatio2s3t, fRatio2s2t, fRatio2s1t, fRatio2s3p5t, fRatio2s2p5t, fRatio2s1p5t, fRatio2s0p5t;
  std::vector<Double_t> fRatio3sfv, fRatio3s3t, fRatio3s2t, fRatio3s1t, fRatio3s3p5t, fRatio3s2p5t, fRatio3s1p5t, fRatio3s0p5t;
  std::vector<Double_t> fHitEffN, fHitEffK;
  Double_t fSpecActivCV, fSpecActivError, fActivCV, fActivError;
  TString fActivID;
  
  TString fGroup;
  TString fFileName;

  ExcelTableValues(TString pdf = "", TString component = "", TString isotope = "", TString mc_id = "") {
    fPdf = pdf;
    fComponent = component;
    fIsotope = isotope;
    fMC_ID = mc_id;
  }

  int FindSuffixIdx(TString suffix){
    std::vector<TString>::iterator pos = std::find(fSuffixes.begin(),fSuffixes.end(),suffix);
    if(pos == fSuffixes.end())
      return -1;
    int si = pos - fSuffixes.begin();
    return si;
  }
  
  void AddSuffix(TString suffix){
    fSuffixes.push_back(suffix);

    fCountsCV.push_back(0.);
    fCountsError.push_back(0.);
    fCountsUL.push_back(0.);

    fRatioFWHMfv.push_back(0.); fRatioFWHM3t.push_back(0.); fRatioFWHM2t.push_back(0.); fRatioFWHM1t.push_back(0.); fRatioFWHM3p5t.push_back(0.); fRatioFWHM2p5t.push_back(0.); fRatioFWHM1p5t.push_back(0.); fRatioFWHM0p5t.push_back(0.);
    fRatio1sfv.push_back(0.); fRatio1s3t.push_back(0.); fRatio1s2t.push_back(0.); fRatio1s1t.push_back(0.); fRatio1s3p5t.push_back(0.); fRatio1s2p5t.push_back(0.); fRatio1s1p5t.push_back(0.); fRatio1s0p5t.push_back(0.);
    fRatio2sfv.push_back(0.); fRatio2s3t.push_back(0.); fRatio2s2t.push_back(0.); fRatio2s1t.push_back(0.); fRatio2s3p5t.push_back(0.); fRatio2s2p5t.push_back(0.); fRatio2s1p5t.push_back(0.); fRatio2s0p5t.push_back(0.);
    fRatio3sfv.push_back(0.); fRatio3s3t.push_back(0.); fRatio3s2t.push_back(0.); fRatio3s1t.push_back(0.); fRatio3s3p5t.push_back(0.); fRatio3s2p5t.push_back(0.); fRatio3s1p5t.push_back(0.); fRatio3s0p5t.push_back(0.);

    fHitEffN.push_back(0.);
    fHitEffK.push_back(0.);
  }

  int FindSuffixIdxAdd(TString suffix){
    if(FindSuffixIdx(suffix) < 0)
      AddSuffix(suffix);
    return FindSuffixIdx(suffix);
  }
  
  void SetHalflife(Double_t hl){
    fHalflife = hl;
  }
  
  void SetExpectedCounts(Double_t counts, Double_t error, Double_t ul, TString suffix){
    int si = FindSuffixIdxAdd(suffix);
    
    fCountsCV[si] = counts;
    fCountsError[si] = error;
    fCountsUL[si] = ul;
  }
  
  void SetHitEfficiency(Double_t n, Double_t k, TString suffix){
    int si = FindSuffixIdxAdd(suffix);
    
    fHitEffN[si] = n;
    fHitEffK[si] = k;
  }

  void SetHitEfficiencyROI(int roi, Double_t fwhmfv, Double_t fwhm3t, Double_t fwhm2t, Double_t fwhm1t, Double_t fwhm3p5t, Double_t fwhm2p5t, Double_t fwhm1p5t, Double_t fwhm0p5t, TString suffix){
    int si = FindSuffixIdxAdd(suffix);

    if(roi == 0)
    {
      fRatioFWHMfv[si] = fwhmfv; fRatioFWHM3t[si] = fwhm3t; fRatioFWHM2t[si] = fwhm2t; fRatioFWHM1t[si] = fwhm1t; fRatioFWHM3p5t[si] = fwhm3p5t; fRatioFWHM2p5t[si] = fwhm2p5t; fRatioFWHM1p5t[si] = fwhm1p5t; fRatioFWHM0p5t[si] = fwhm0p5t;
    }
    else if(roi == 1)
    {
      fRatio1sfv[si] = fwhmfv; fRatio1s3t[si] = fwhm3t; fRatio1s2t[si] = fwhm2t; fRatio1s1t[si] = fwhm1t; fRatio1s3p5t[si] = fwhm3p5t; fRatio1s2p5t[si] = fwhm2p5t; fRatio1s1p5t[si] = fwhm1p5t; fRatio1s0p5t[si] = fwhm0p5t;
    }
    else if(roi == 2)
    {
      fRatio2sfv[si] = fwhmfv; fRatio2s3t[si] = fwhm3t; fRatio2s2t[si] = fwhm2t; fRatio2s1t[si] = fwhm1t; fRatio2s3p5t[si] = fwhm3p5t; fRatio2s2p5t[si] = fwhm2p5t; fRatio2s1p5t[si] = fwhm1p5t; fRatio2s0p5t[si] = fwhm0p5t;
    }
    else if(roi == 3)
    {
      fRatio3sfv[si] = fwhmfv; fRatio3s3t[si] = fwhm3t; fRatio3s2t[si] = fwhm2t; fRatio3s1t[si] = fwhm1t; fRatio3s3p5t[si] = fwhm3p5t; fRatio3s2p5t[si] = fwhm2p5t; fRatio3s1p5t[si] = fwhm1p5t; fRatio3s0p5t[si] = fwhm0p5t;
    }
  }

  void SetActivity(Double_t specActiv, Double_t specActivError, Double_t activ, Double_t activError, TString activID){
    fSpecActivCV = specActiv;
    fSpecActivError = specActivError;
    fActivCV = activ;
    fActivError = activError;
    fActivID = activID;
  }    

  void SetGroup(TString group){
    fGroup = group;
  }

  void SetFileName(TString filename){
    fFileName = filename;
  }
  
  void Copy(ExcelTableValues& other){
    fPdf = other.fPdf;
    fComponent = other.fComponent;
    fIsotope = other.fIsotope;

    fHalflife = other.fHalflife;
    
    fSuffixes.clear();
    
    fCountsCV.clear();
    fCountsError.clear();
    fCountsUL.clear();

    fRatioFWHMfv.clear(); fRatioFWHM3t.clear(); fRatioFWHM2t.clear(); fRatioFWHM1t.clear(); fRatioFWHM3p5t.clear(); fRatioFWHM2p5t.clear(); fRatioFWHM1p5t.clear(); fRatioFWHM0p5t.clear();
    fRatio1sfv.clear(); fRatio1s3t.clear(); fRatio1s2t.clear(); fRatio1s1t.clear(); fRatio1s3p5t.clear(); fRatio1s2p5t.clear(); fRatio1s1p5t.clear(); fRatio1s0p5t.clear();
    fRatio2sfv.clear(); fRatio2s3t.clear(); fRatio2s2t.clear(); fRatio2s1t.clear(); fRatio2s3p5t.clear(); fRatio2s2p5t.clear(); fRatio2s1p5t.clear(); fRatio2s0p5t.clear();
    fRatio3sfv.clear(); fRatio3s3t.clear(); fRatio3s2t.clear(); fRatio3s1t.clear(); fRatio3s3p5t.clear(); fRatio3s2p5t.clear(); fRatio3s1p5t.clear(); fRatio3s0p5t.clear();

    fHitEffN.clear();
    fHitEffK.clear();

    for(size_t i = 0; i < other.fSuffixes.size(); i++)
    {
      fSuffixes.push_back(other.fSuffixes.at(i));

      if(i < other.fCountsCV.size())
      {
        fCountsCV.push_back(other.fCountsCV.at(i));
        fCountsError.push_back(other.fCountsError.at(i));
        fCountsUL.push_back(other.fCountsUL.at(i));
      }
      
      if(i < other.fHitEffN.size())
      {
        fHitEffN.push_back(other.fHitEffN.at(i));
        fHitEffK.push_back(other.fHitEffK.at(i));
        fRatioFWHMfv.push_back(other.fRatioFWHMfv.at(i)); fRatioFWHM3t.push_back(other.fRatioFWHM3t.at(i)); fRatioFWHM2t.push_back(other.fRatioFWHM2t.at(i)); fRatioFWHM1t.push_back(other.fRatioFWHM1t.at(i)); fRatioFWHM3p5t.push_back(other.fRatioFWHM3p5t.at(i)); fRatioFWHM2p5t.push_back(other.fRatioFWHM2p5t.at(i)); fRatioFWHM1p5t.push_back(other.fRatioFWHM1p5t.at(i)); fRatioFWHM0p5t.push_back(other.fRatioFWHM0p5t.at(i));
        fRatio1sfv.push_back(other.fRatio1sfv.at(i)); fRatio1s3t.push_back(other.fRatio1s3t.at(i)); fRatio1s2t.push_back(other.fRatio1s2t.at(i)); fRatio1s1t.push_back(other.fRatio1s1t.at(i)); fRatio1s3p5t.push_back(other.fRatio1s3p5t.at(i)); fRatio1s2p5t.push_back(other.fRatio1s2p5t.at(i)); fRatio1s1p5t.push_back(other.fRatio1s1p5t.at(i)); fRatio1s0p5t.push_back(other.fRatio1s0p5t.at(i));
        fRatio2sfv.push_back(other.fRatio2sfv.at(i)); fRatio2s3t.push_back(other.fRatio2s3t.at(i)); fRatio2s2t.push_back(other.fRatio2s2t.at(i)); fRatio2s1t.push_back(other.fRatio2s1t.at(i)); fRatio2s3p5t.push_back(other.fRatio2s3p5t.at(i)); fRatio2s2p5t.push_back(other.fRatio2s2p5t.at(i)); fRatio2s1p5t.push_back(other.fRatio2s1p5t.at(i)); fRatio2s0p5t.push_back(other.fRatio2s0p5t.at(i));
        fRatio3sfv.push_back(other.fRatio3sfv.at(i)); fRatio3s3t.push_back(other.fRatio3s3t.at(i)); fRatio3s2t.push_back(other.fRatio3s2t.at(i)); fRatio3s1t.push_back(other.fRatio3s1t.at(i)); fRatio3s3p5t.push_back(other.fRatio3s3p5t.at(i)); fRatio3s2p5t.push_back(other.fRatio3s2p5t.at(i)); fRatio3s1p5t.push_back(other.fRatio3s1p5t.at(i)); fRatio3s0p5t.push_back(other.fRatio3s0p5t.at(i));
      }
    }

    fSpecActivCV = other.fSpecActivCV;
    fSpecActivError = other.fSpecActivError;
    fActivCV = other.fActivCV;
    fActivError = other.fActivError;
    fActivID = other.fActivID;

    fGroup = other.fGroup;
    fFileName = other.fFileName;
  }

  Double_t GetMeanCounts(TString suffix, bool eval = false){

    int si = FindSuffixIdx(suffix);
    if(si < 0)
      return -1;
          
    if(not eval and si < (int) fCountsCV.size())
      return fCountsCV[si];

    if(si > (int) fHitEffK.size())
      return -1;
    
    double time = 1; // years
    if(time/fHalflife > 0.01)
    {
      double lhl = TMath::Log2(fHalflife);
      time = lhl * (1 - exp(-time/lhl));
    }

    double eff = fHitEffK[si]/fHitEffN[si];
    double activ = fActivCV;

    double counts = time * eff * activ * 31556736; // seconds per year conversion

    return counts;
  }

  Double_t GetEfficiency(){

    Double_t totEff = std::accumulate(fHitEffK.begin(),fHitEffK.end(),0.);
    Double_t totN = *std::max_element(fHitEffN.begin(),fHitEffN.end());

    return totN > 0 ? totEff*1./totN : -1;
  }

  Double_t GetEfficiencyError(){

    Double_t totEff = std::accumulate(fHitEffK.begin(),fHitEffK.end(),0.);
    Double_t totN = *std::max_element(fHitEffN.begin(),fHitEffN.end());

    return totN > 0 ? sqrt(totEff)/totN : -1;
  }  

  TString GetMC_ID(){
     return fMC_ID;
  }

  void Print(){
    std::cout << Form("* PDF : %s , Group: %s ", fPdf.Data(), fGroup.Data()) << std::endl;
    std::cout << Form("** Component: %s , Isotope: %s , Halflife: %g", fComponent.Data(), fIsotope.Data(), fHalflife) << std::endl;
    for(size_t i = 0; i < fSuffixes.size(); i++)
    {
      if(i < fCountsCV.size())
        std::cout << Form("** %s Counts: CV = %g , Error = %g , UL = %g ", fSuffixes[i].Data(), fCountsCV[i], fCountsError[i], fCountsUL[i]) << std::endl;
      if(i < fHitEffN.size())
        std::cout << Form("** %s Hit efficiency: N = %g , k = %g , FWHM FV ratio: %g , FWHM 3t ratio: %g , 2-sigma 2t ratio: %g , FWHM 1t ratio: %g", fSuffixes[i].Data(), fHitEffN[i], fHitEffK[i], fRatioFWHMfv[i], fRatioFWHM3t[i], fRatio2s2t[i], fRatioFWHM1t[i]) << std::endl;
      //std::cout << Form("** %s Internal Counts = %g", fSuffixes[i].Data(), GetMeanCounts(fSuffixes[i],true)) << std::endl;
    }
    std::cout << Form("** Activity: Specific CV = %g , Specific Error = %g , Full CV = %g , Full Error = %g , ID = %s", fSpecActivCV, fSpecActivError, fActivCV, fActivError, fActivID.Data()) << std::endl;
    std::cout << Form("** File: %s", fFileName.Data()) << std::endl;
  }
      
  
} ExcelTableValues;


class nEXOSensitivity : public TObject
{
public:

  enum ExpectCount {kPosCV = 0, kUL, kPosUL, kRdmCV};

  nEXOSensitivity(int seed = 0, const char* treeFileName = 0);
  virtual ~nEXOSensitivity();

  void SetSeed(int seed);
  void SetBinning(int nBinsX, double xMin, double xMax, int nBinsY, double yMin, double yMax);
    
  void SetBaTag(bool useBaTag){fBaTag = useBaTag; if(useBaTag) fWithStandoff = false;};
  void TurnGroupOff(TString groupName){fTurnedOffGroups.insert(groupName);}
  void ScaleMeanBackgrounds(Double_t scale, Bool_t bb2n){fScaleBkgds = scale; fScale2nu = bb2n;}
  
  void AddUserMeanCounts(TString group, Double_t value){
    fUserMeanCounts.insert(std::make_pair(group.Data(),value));
  };

  void SetNoFracFlag( bool noFracFlag ){ fNoFracFlag = noFracFlag; }
  
  void MakeFittingHistogramFile();
  void BuildWorkspace(Double_t yrs = 5.0, Double_t signalCounts = 0.0);
  void AddMagicNumber(double signal, double magicN);
  void GenAndFitData(Int_t nRuns = 1, Double_t yrs = 5., Double_t signalCounts = 0., Int_t rdmRate = 0);

  void LoadComponentHistograms();
  void MakeGroupHistograms();

  Double_t EvalCounts(Double_t hitEfficiency, Double_t activity, Double_t time, Double_t halflife);
  Double_t GenTruncGaus(Double_t mean, Double_t sigma, Double_t unif, Int_t low = 0.);

  TGraph * GetLambdaVsHypGraph() { return g_nllratio_vs_hyp; }


  TString fResultFileName;

  ExpectCount fExpectCountMethod;
  Bool_t fRunTruthValFit;
  Int_t fVerboseLevel;
  Double_t fSSFracImprovement;
  Double_t fRn222RateCorrection;
      
protected:
    
  TString fTreeFileName;
  TTree* fExcelTree;
  
  RooWorkspace* fWsp;
  TString fWspFileName;
  bool fWriteWsp;
   
  
  std::vector<std::unique_ptr<RooDataHist> > hdata_ss;
  std::vector<std::unique_ptr<RooDataHist> > hdata_ms;
  std::vector<std::unique_ptr<RooDataHist> > hdata1D_ss;
  std::vector<std::unique_ptr<RooDataHist> > hdata1D_ms;
  std::vector<std::unique_ptr<RooRealVar> > num_tomakefit; //number of events for component (floating)
  std::vector<std::unique_ptr<RooRealVar> > frac_tomakefit; //ss fraction for component (floating)
  std::vector<std::unique_ptr<RooFormulaVar> > coef_tomakefit; //coefficient in sum pdf for component
  
  
  std::map<TString, TH1*> fComponentHistos;
  std::map<TString, TH1*> fGroupHistos;
  std::map<TString, double> fNormHistoBins;

  std::map<TString, Double_t> fUserMeanCounts;

  Bool_t fBaTag;
  std::set<TString> fTurnedOffGroups;
  Double_t fScaleBkgds;
  Bool_t fScale2nu;
  
  Bool_t fRandomizeMeanCounts;
  Bool_t fWithStandoff;

  void LoadExcelTree(const char* filename,const char* treename = "ExcelTableValues");
  void ReadExcelTree();
  void SetUserMeanCounts(TString groupName, Double_t value);
  void SetAllGroupMeanCounts(Double_t value, TString except = "LXeBb2n");

  // MakeFittingHistogramFile
  TH2D* AdjustedBinHist(TH2D& inHist);
  TH2D* MakeCombinedHisto(TString histName, Int_t nComp, TString* compNames, Double_t* fractions, TString isotope, Bool_t isSS);

  // BuildWorkspace
  void BuildGenHistos(TString* pdfNames, const int nPdfs, Double_t* meanPerYear_ss, Double_t* meanPerYear_ms, Double_t yrs, Double_t signalCounts);
  void BuildFitPdf(TString* pdfNames, const int nFitPdfs, Bool_t isSS);
  void BuildFitPdfs(TString* pdfNames, const int nFitPdfs);

  // GenAndFitData
  void GetFitPdfNames(RooAddPdf* fitPdf);
  RooGaussian* GetRn222Constraint(Double_t rateError, Bool_t randomize);
  RooGaussian* GetEfficiencyConstraint(Double_t effError, Bool_t randomize);
  RooMultiVarGaussian* GetFracConstraint(Double_t fracError, Bool_t randomize);
  RooAbsData* GenerateData(RooAbsPdf* genPdf, RooArgSet obs, Bool_t isBinned);
  std::map<double,double> fMagic_numbers;
  // Member variables
  Int_t fSeed;
  TRandom3 fRandom;
  TGraph * g_nllratio_vs_hyp;
  
  //std::vector<TString> fAllComponents;
  std::vector<TString> fComponentNames;
  std::map<TString, std::vector<TString> > fGroups;
  std::map<TString, Double_t > fGroupMeanCounts;  
  std::map<TString, std::vector<Double_t> > fGroupFractions;
  std::map<TString, std::pair<Double_t,Double_t> > fGroupSeparation;
  
  std::map<TString, Double_t > fGroupMeanRatioFVfwhm;  
  std::map<TString, Double_t > fGroupMeanRatio3tfwhm;  
  std::map<TString, Double_t > fGroupMeanRatio2tfwhm;  
  std::map<TString, Double_t > fGroupMeanRatio1tfwhm;  
  std::map<TString, Double_t > fGroupMeanRatio3p5tfwhm;  
  std::map<TString, Double_t > fGroupMeanRatio2p5tfwhm;  
  std::map<TString, Double_t > fGroupMeanRatio1p5tfwhm;  
  std::map<TString, Double_t > fGroupMeanRatio0p5tfwhm;  

  std::map<TString, Double_t > fGroupMeanRatioFV1sig;  
  std::map<TString, Double_t > fGroupMeanRatio3t1sig;  
  std::map<TString, Double_t > fGroupMeanRatio2t1sig;  
  std::map<TString, Double_t > fGroupMeanRatio1t1sig;  
  std::map<TString, Double_t > fGroupMeanRatio3p5t1sig;  
  std::map<TString, Double_t > fGroupMeanRatio2p5t1sig;  
  std::map<TString, Double_t > fGroupMeanRatio1p5t1sig;  
  std::map<TString, Double_t > fGroupMeanRatio0p5t1sig;  

  std::map<TString, Double_t > fGroupMeanRatioFV2sig;  
  std::map<TString, Double_t > fGroupMeanRatio3t2sig;  
  std::map<TString, Double_t > fGroupMeanRatio2t2sig;  
  std::map<TString, Double_t > fGroupMeanRatio1t2sig;  
  std::map<TString, Double_t > fGroupMeanRatio3p5t2sig;  
  std::map<TString, Double_t > fGroupMeanRatio2p5t2sig;  
  std::map<TString, Double_t > fGroupMeanRatio1p5t2sig;  
  std::map<TString, Double_t > fGroupMeanRatio0p5t2sig;  

  std::map<TString, Double_t > fGroupMeanRatioFV3sig;  
  std::map<TString, Double_t > fGroupMeanRatio3t3sig;  
  std::map<TString, Double_t > fGroupMeanRatio2t3sig;  
  std::map<TString, Double_t > fGroupMeanRatio1t3sig;  
  std::map<TString, Double_t > fGroupMeanRatio3p5t3sig;  
  std::map<TString, Double_t > fGroupMeanRatio2p5t3sig;  
  std::map<TString, Double_t > fGroupMeanRatio1p5t3sig;  
  std::map<TString, Double_t > fGroupMeanRatio0p5t3sig;  

  int fNbinsX;
  double fXmin;
  double fXmax;
  int fNbinsY;
  double fYmin;
  double fYmax;
  std::map<int,int> fbinAdjustMap;


  Double_t fBkgdTotal, fBkgdFwhmFV, fBkgdFwhm3t, fBkgdFwhm2t, fBkgdFwhm1t, fBkgdFwhm3p5t, fBkgdFwhm2p5t, fBkgdFwhm1p5t, fBkgdFwhm0p5t;
  Double_t fBkgd1sigmaFV, fBkgd1sigma3t, fBkgd1sigma2t, fBkgd1sigma1t, fBkgd1sigma3p5t, fBkgd1sigma2p5t, fBkgd1sigma1p5t, fBkgd1sigma0p5t;
  Double_t fBkgd2sigmaFV, fBkgd2sigma3t, fBkgd2sigma2t, fBkgd2sigma1t, fBkgd2sigma3p5t, fBkgd2sigma2p5t, fBkgd2sigma1p5t, fBkgd2sigma0p5t;
  Double_t fBkgd3sigmaFV, fBkgd3sigma3t, fBkgd3sigma2t, fBkgd3sigma1t, fBkgd3sigma3p5t, fBkgd3sigma2p5t, fBkgd3sigma1p5t, fBkgd3sigma0p5t;

  Bool_t fWithEff;
  Int_t fInterpOrder;
  Double_t fMeanSignalEff;
  Double_t fSignalEffError;
  Double_t fFracError;
  Double_t fRnRateError;
  TString fSignalName;
  Double_t fFidVol;
  Bool_t fNoFracFlag;

  Int_t fNcpu; 							//number of cpus to parallelize over
  std::vector<TString> fFitPdfNames;	//a vector list of the fitting pdf names
  int fNFitPdfs;							//number of fitting pdfs
  double fErrorLevel;		//minuit error level, to get the 90% UL
  int fPrintLevel;				//minuit print level

  ClassDef(nEXOSensitivity,1)
};

#endif


