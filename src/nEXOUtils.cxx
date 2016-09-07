//////////////////////////////////////////////////////////////////////////////////
//
//
// nEXOUtils is a namespace with a collection of functions useful for sensitivity
//
//
//////////////////////////////////////////////////////////////////////////////////

#include "nEXOUtils.hh"

Double_t nEXOUtils::GetNumXeAtoms(Double_t xeMass, Double_t enrichment)
{
  // return the number of atoms in a given xenon mass and enrichment
  
  return xeMass * enrichment * AVOGADRO_NUMBER / XE136_MOLAR_MASS;
}

Double_t nEXOUtils::GetHalfLife(Double_t counts, Double_t yrs, Double_t xeMass)
{
  // return half-life for given counts, years and mass

  return GetNumXeAtoms(xeMass) * yrs * log(2.) / counts;
}

Double_t nEXOUtils::GetSensHalfLife(TString inFileName, Double_t yrs, Double_t xeMass, TString signalName, double eff)
{
  // return the half-life estimated from the upper limit counts of the tree in given file

  double counts = GetBackgroundCounts(inFileName,Form("num_%s+num_%s_eHi",signalName.Data(),signalName.Data()));
  counts /= eff;
  double hl = GetHalfLife(counts,yrs,xeMass);
  
  return hl;
}
  
Double_t nEXOUtils::GetBackgroundCounts(TString inFileName, TString varName, TString treeName, bool useMean)
{
  // return the median counts from the tree named "treeName" using variable named "varName"
  // current var names: bkg_fwhm_fv, bkg_fwhm_3t, bkg_fwhm_1t, roi_bkg, roi_bkg_3t
  
  TChain tree(treeName.Data());
  tree.Add(inFileName.Data());

  tree.SetEstimate(tree.GetEntries()+1);
  tree.Draw(varName.Data(),"covQual_sig == 3 && stat_sig == 0","goff");

  Double_t median = GetBackgroundCounts(tree.GetSelectedRows(),tree.GetV1(),useMean);

  return median;  
}

Double_t nEXOUtils::GetBackgroundCounts(Long64_t length, Double_t* counts,bool useMean)
{
  if(useMean)
    return TMath::Mean(length,counts);

  return TMath::Median(length,counts);
}

Double_t nEXOUtils::GetDiscoveryCutValue(TH1D* histo, Double_t gausSignif)
{
  if (!histo) {std::cout<<"Histo is empty!"<<std::endl; return -1.;}
  Double_t alpha = ROOT::Math::normal_cdf_c(gausSignif);
  
  Int_t nBins = histo->GetNbinsX();
  Double_t integral = histo->Integral(1, nBins);
  Double_t sum = 0.;
  Int_t index = -1;
  while (sum < alpha) {
    index++;
    sum += histo->GetBinContent(nBins-index)/integral;
  }
  return histo->GetXaxis()->GetBinLowEdge(nBins-index);
}


bool nEXOUtils::BuildNLLRatioHistogram(TH1D& histo, const char* filename, const char* treename, const char* cut)
{
  TChain chain(treename);
  chain.Add(filename);

  chain.SetEstimate(chain.GetEntries()+1);
  chain.Draw("nll_ratio:1",cut,"goff");

  if(chain.GetSelectedRows() <= 0)
    return false;
  
  histo.FillN(chain.GetSelectedRows(),chain.GetV1(),chain.GetV2());
  histo.Scale(1./histo.Integral());
  
  histo.GetXaxis()->SetTitle("NLL Ratio");
  
  return true;
}

Double_t nEXOUtils::GetDiscoveryCounts(Double_t prob, Double_t sigma, Long64_t length, Double_t* counts, Double_t yrs, const char* filename, double max, double min, int bins, const char* treename, const char* cut)
{
  std::map<TString, TH1D> nllRatioHistos;
  
  for(int i = 0; i < length; i++)
  {
    double count = counts[i];
    TString file = Form(filename,yrs,count);
    TH1D nllRatioHisto(Form("nllRatioHisto_%0.1f",count),"",bins,min,max);

    BuildNLLRatioHistogram(nllRatioHisto,file,treename,cut);

    nllRatioHistos.insert(std::make_pair<TString,TH1D>(Form("%0.1f",count),TH1D(nllRatioHisto)));
  }

  Double_t cutVal = GetDiscoveryCutValue(&nllRatioHistos.at("0.0"),sigma);

  TGraph grCounts;
  for(int i = 0; i < length; i++)
  {
    double count = counts[i];
    TString countName = Form("%0.1f",count);
    if(countName == "0.0")
      continue;
    TH1D* histo = &nllRatioHistos.at(countName);
    Int_t cutBin = histo->GetXaxis()->FindBin(cutVal);

    double y = histo->Integral(cutBin, bins); // histo should be normalized
    //std::cout << "count " << count << " val " << y << std::endl;
    grCounts.SetPoint(grCounts.GetN(), y, count);
  }

  TFile rout(Form("disc_graph_%0.1f_yrs_%0.1f_prob_%0.1f_sigma.root",yrs,prob,sigma),"recreate");
  grCounts.Write("g");
  rout.Close();

  Double_t result = grCounts.Eval(prob);
  std::cout << "Discovery counts = " << result << std::endl;
  
  return result;
}

Double_t nEXOUtils::GetDiscHalfLife(Double_t prob, Double_t sigma, Long64_t length, Double_t* counts, Double_t yrs, const char* filename, double mass, double max, double min, int bins, const char* treename, const char* cut)
{
  // return the half-life estimated from the upper limit counts of the tree in given file

  double count = GetDiscoveryCounts(prob,sigma,length,counts,yrs,filename,max,min,bins,treename,cut);
  double hl = GetHalfLife(count,yrs,mass);
  
  return hl;
}

Double_t* nEXOUtils::ReadSensFromFiles(size_t n, double* yrs, const char* files, double mass)
{
  std::cout << "Reading sensitivity ... \n";
  if(n == 0 or not yrs or not files)
    return 0;

  Double_t* sens = new Double_t[n];
    
  for(size_t i = 0; i < n; i++)
  {
    TString fileName(Form(files,yrs[i]));
    sens[i] = GetSensHalfLife(fileName.Data(),yrs[i],mass);
    std::cout << "HL " << yrs[i] << " = " << sens[i] << std::endl;
  }
  
  return sens;    
}

Double_t* nEXOUtils::ReadDiscFromFiles(size_t ny, double* yrs, size_t nc, double* cts, const char* files, double mass, Double_t prob, Double_t sigma, double max, double min, int bins, const char* treename, const char* cut)
{
  std::cout << "Reading discovery potential ... \n";
  if(ny == 0 or not yrs or nc == 0 or not cts or not files)
    return 0;

  Double_t* disc = new Double_t[ny];
    
  for(size_t i = 0; i < ny; i++)
  {
    disc[i] = GetDiscHalfLife(prob,sigma,nc,cts,yrs[i],files,mass,max,min,bins,treename,cut);
    std::cout << "HL " << yrs[i] << " = " << disc[i] << std::endl;
  }
  
  return disc;    
}

Double_t nEXOUtils::EvalCountSensHalfLife(TString inFileName, TString varName, Double_t yrs, Double_t mass)
{
  Double_t mean = GetBackgroundCounts(inFileName,varName,"tree",false);
  mean *= yrs;
  TFeldmanCousins fc(0.9);
  fc.SetMuMax(mean < 1 ? 50 : 50*mean);
  fc.SetMuMin(0.0);
  fc.SetMuStep((fc.GetMuMax()-fc.GetMuMin())/50000.);
  
  Double_t avg = 0.;
  Int_t obs = 0;
  Double_t diff = 0.;
  do
  {
    Double_t prev = avg;
    Double_t fcul = fc.CalculateUpperLimit(obs,mean);
    Double_t pmo = TMath::PoissonI(obs,mean);
    avg += pmo * fcul;
    //avg += TMath::PoissonI(obs,mean) * fc.CalculateUpperLimit(obs,mean);
    diff = avg - prev;
    //std::cout << Form("mean = %g, obs = %d, pmo = %g, fcul = %g, prev = %g, avg %g, diff = %g",mean,obs,pmo,fcul,prev,avg,diff) << std::endl;
    obs++;
  }while(diff > 0.01 || obs < mean);
  
  return GetHalfLife(avg,yrs,mass);
}

Double_t* nEXOUtils::EvalCountSensFromFiles(size_t n, double* yrs, const char* files, TString varName, double mass, double eff)
{
  std::cout << "Evaluating counting sensitivity ... \n";
  if(n == 0 or not yrs or not files)
    return 0;

  Double_t* sens = new Double_t[n];
    
  for(size_t i = 0; i < n; i++)
  {
    TString fileName(Form(files,yrs[i]));
    sens[i] = eff * EvalCountSensHalfLife(fileName.Data(),varName,yrs[i],mass);
    std::cout << "Counting HL " << yrs[i] << " = " << sens[i] << " in mass = " << mass << std::endl;
  }
  
  return sens;    
   
}

Double_t nEXOUtils::GetCriticalCounts(Double_t alpha, Double_t nBkgd)
{
  //Takes in the required significance and background counts and returns the critical number of counts
  Double_t prob = 1.;
  Int_t counts = 0;
  while (prob > alpha) {
    prob = ROOT::Math::poisson_cdf_c(counts, nBkgd);
    counts++;
  }
  return (Double_t)counts;
}

Double_t nEXOUtils::GetLeastDetectableSignal(Double_t prob, Double_t nBkgd, Double_t nCrit)
{
  //Takes in the required probability, background counts, and critical counts and returns the least detectable signal counts
  TF1 f("func", Form("ROOT::Math::poisson_cdf_c(%i, x)-%0.2f", (Int_t)nCrit, prob), 0., 1000.);
  ROOT::Math::WrappedTF1 wf1(f);
  
  ROOT::Math::BrentRootFinder brf;
  brf.SetNpx(50000);
  brf.SetFunction(wf1, 0., 1000.);
  brf.Solve();
  
  return brf.Root()-nBkgd;
}

Double_t nEXOUtils::EvalCountDiscHalfLife(TString inFileName, TString varName, Double_t yrs, Double_t mass)
{
  Double_t mean = GetBackgroundCounts(inFileName,varName,"tree");
  mean *= yrs;

  Double_t significance = 3.; // sigmas
  Double_t alpha = ROOT::Math::normal_cdf_c(significance, 1.);
  Double_t cl = 0.5;
  
  Double_t critical = GetCriticalCounts(alpha,mean);
  Double_t counts = GetLeastDetectableSignal(cl, mean, critical);
    
  return GetHalfLife(counts,yrs,mass);
}

Double_t* nEXOUtils::EvalCountDiscFromFiles(size_t n, double* yrs, const char* files, TString varName, double mass, double eff)
{
  std::cout << "Evaluating counting discovery ... \n";
  if(n == 0 or not yrs or not files)
    return 0;

  Double_t* disc = new Double_t[n];
    
  for(size_t i = 0; i < n; i++)
  {
    TString fileName(Form(files,yrs[i]));
    disc[i] = eff * EvalCountDiscHalfLife(fileName.Data(),varName,yrs[i],mass);
    std::cout << "Counting discovery HL " << yrs[i] << " = " << disc[i] << " in mass = " << mass << std::endl;
  }
  
  return disc;    
}
