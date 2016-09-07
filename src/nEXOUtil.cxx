#include "nEXOUtil.hh"

ClassImp(nEXOUtil);

//-----------------------------------------------------------
nEXOUtil::nEXOUtil() {
	Double_t matElem[5] = {0.31, 2.21, 0.83, 1.78, 0.49};
	TString	matElemName[5] = {"gcm", "qrpa2", "rqrpa", "nsm", "ibm2"};
	TString	matElemTitle[5] = {"GCM", "QRPA-2", "R-QRPA", "NSM", "IBM-2"};
	TString	matElemAxisTitle[5] = {"m_{#beta#beta} (meV) GCM", "m_{#beta#beta} (meV) QRPA-2", "m_{#beta#beta} (meV) R-QRPA", "m_{#beta#beta} (meV) NSM", "m_{#beta#beta} (meV) IBM-2"};
	for (Int_t i=0; i<5; i++) {
		fMatElem[i] = matElem[i];
		fMatElemName[i] = matElemName[i];
		fMatElemTitle[i] = matElemTitle[i];
		fMatElemAxisTitle[i] = matElemAxisTitle[i];
	}
	fXaxisTitle = "Exposure (yrs)";
	fYaxisTitle = "^{136}Xe T_{1/2}";

	fFidVolXY = 628.;
	fFidVolZZ = 631.;
}

//-----------------------------------------------------------
nEXOUtil::~nEXOUtil() {}


//-----------------------------------------------------------
double nEXOUtil::getMedian(TH1D* h1, double& lo, double& hi)
{
  int n = h1->GetXaxis()->GetNbins(); 
  double integral = 0.;
  int index = 0;
  double sum = h1->Integral(1, n);
  while (integral < 0.5) {
    index++;
    integral += h1->GetBinContent(index)/sum;
  }
  double loEdge = h1->GetXaxis()->GetBinLowEdge(index);
  double binWidth = h1->GetXaxis()->GetBinWidth(1);
  double frac = integral - 0.5;
  
  double median = loEdge+frac*binWidth;
  double target = (ROOT::Math::normal_cdf(1)-ROOT::Math::normal_cdf(-1))/2.;
  integral = h1->GetBinContent(index)/sum*frac;
  int index_lo = index;
  while (integral < target/2.) {
    index_lo--;
    integral += h1->GetBinContent(index_lo)/sum;
    std::cout << integral << std::endl;
  }
  double loEdge_lo = h1->GetXaxis()->GetBinLowEdge(index_lo);
  double frac_lo = TMath::Abs(integral - target);
  lo = loEdge_lo;
  
  integral = h1->GetBinContent(index)/sum*frac;
  int index_hi = index;
  while (integral < target/2.) {
    index_hi++;
    integral += h1->GetBinContent(index_hi)/sum;
  }
  double loEdge_hi = h1->GetXaxis()->GetBinUpEdge(index_hi);
  double frac_hi = TMath::Abs(integral - target);
  hi = loEdge_hi;
  
  return median;
}

//-----------------------------------------------------------
void nEXOUtil::PlotBb0nUL(TString inFileName, TString outFileName, Int_t withBb0n)
{
  TFile fIn(inFileName.Data());
  TTree* tree = (TTree*)fIn.Get("tree");
  
  TFile fOut(outFileName.Data(), "recreate");
  TH1D* h_num = new TH1D("h_num", "Best Fit #beta#beta0#nu Counts", 200, 0., 200.);
  TH1D* h_ehi = new TH1D("h_ehi", "", 100, 0., 100.);		
  TH1D* h_elo = new TH1D("h_elo", "", 100, 0., 100.);
  TH1D* h_ul	= new TH1D("h_up","90% UL on #beta#beta0#nu Counts", 1000, 0., 50.);		
  
  h_num->GetXaxis()->SetTitle("#beta#beta0#nu Counts (Best fit)");
  h_ul->GetXaxis()->SetTitle("Upper Limit on #beta#beta0#nu Counts (90% CL)");
  
  //set the addresses of the wanted variables
  Double_t num_bb0n = 0.;
  Double_t num_bb0n_eLo = 0.;
  Double_t num_bb0n_eHi = 0.;
  Double_t num_bb0n_ul = 0.;
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("num_LXeBb0n", 1);
  tree->SetBranchStatus("num_LXeBb0n_eLo", 1);
  tree->SetBranchStatus("num_LXeBb0n_eHi", 1);
  //tree->SetBranchStatus("num_bb0n_ul", &num_bb0n_ul);
  tree->SetBranchAddress("num_LXeBb0n", &num_bb0n);
  tree->SetBranchAddress("num_LXeBb0n_eLo", &num_bb0n_eLo);
  tree->SetBranchAddress("num_LXeBb0n_eHi", &num_bb0n_eHi);
  //tree->SetBranchAddress("num_bb0n_ul", &num_bb0n_ul);
  
  //tree->SetBranchAddress("fitres_sig",&fitres);
  
  //Fill the histograms
  for (int i=0; i<tree->GetEntries(); i++) {
    //cout << i << " " << tree->GetEntries() << endl;
    tree->GetEntry(i);
    h_num->Fill(num_bb0n);
    h_elo->Fill(-num_bb0n_eLo);
    h_ehi->Fill(num_bb0n_eHi);
    h_ul->Fill(num_bb0n+num_bb0n_eHi);	
    
  }
  
  //construct 2D color plot of correlation matrix
  //gStyle->SetOptStat(0);
  //gSTyle->SetPalette(1);
  //TH2* hcorr = tree->CorrelationHist();
  
  // Extract correlation matrix as TMatrixDSym
  //const TMatrixDSym& cor = tree->correlationMatrix();
  //Print correlation matrix
  //cout<<"correlation matrix" << endl;
  //cor.Print();
  
  //calculate the mean, median, and 68% interval containing the median
  Double_t loE = 0., hiE = 0.;
  Double_t mean = h_ul->GetMean();
  Double_t median = getMedian(h_ul, loE, hiE);
  std::cout << "Mean counts = " << mean << ", median 90% UL = " << median << std::endl;
  std::cout << loE << " , " << hiE << std::endl;
  
  
  //Find the bins for the 68% interval
  h_ul->Rebin(20);
  Int_t loBin = h_ul->GetXaxis()->FindBin(loE);
  Int_t hiBin = h_ul->GetXaxis()->FindBin(hiE);
  
  //Histo to display the 68% interval around the median
  TH1D* h_interval = new TH1D("h_interval", "", 100, 0., 100.);
  for (int i=loBin; i<=hiBin; i++) {
    h_interval->SetBinContent(i, h_ul->GetBinContent(i));
  }
  h_interval->SetFillColor(kCyan);
  h_interval->SetFillStyle(3001);
  h_interval->SetLineColor(0);
  
  //Find the bins corresponding to the median upper limit, and the upper limit from real data
  Double_t medBinCon = h_ul->GetBinContent(h_ul->GetXaxis()->FindBin(median));
  //Double_t fitBinCon = h_ul->GetBinContent(h_ul->GetXaxis()->FindBin(fitUL));	
  
  //Draw the results
  TCanvas* cc = new TCanvas("cc","");
  cc->Divide(2,1);
  cc->cd(1);
  h_num->DrawCopy();
  cc->cd(2);
  h_ul->DrawCopy();
  h_interval->DrawCopy("same");
  h_ul->DrawCopy("same");
  TLine* medianLine = new TLine(median, 0., median, medBinCon);
  medianLine->SetLineColor(kRed);
  medianLine->SetLineWidth(2);
  //TLine* fitMedianLine = new TLine(fitUL, 0., fitUL, fitBinCon);
  //fitMedianLine->SetLineColor(kBlue);
  //fitMedianLine->SetLineWidth(2);
  
  TLegend* leg = new TLegend(0.475, 0.45, 0.925, 0.7);
  //leg->AddEntry(fitMedianLine, "Upper limit from data, 23.92 counts", "l");
  TString entryTitle = "Median upper limit"; entryTitle += Form(", %0.2f counts", median);
  leg->AddEntry(medianLine, entryTitle, "l");
  leg->AddEntry(h_interval, "68% Interval around Median", "f");
  
  medianLine->Draw();
  //fitMedianLine->Draw();
  leg->Draw();
  gPad->RedrawAxis();
  cc->SaveAs("cc.C");
  
  TCanvas* cc2 = new TCanvas("cc2", "");
  h_ul->DrawCopy();
  h_interval->DrawCopy("same");
  h_ul->DrawCopy("same");
  medianLine->Draw();
  //fitMedianLine->Draw();
  leg->Draw();
  gPad->RedrawAxis();	
  
  fOut.cd();
  cc->Write();
  h_num->Write();
  h_elo->Write();
  h_ehi->Write();
  h_ul->Write();
}

//-----------------------------------------------------------
double nEXOUtil::getHalfLife(TString inFileName, Double_t yrs, Double_t xeMass)
{
  TChain tree("tree");
  tree.AddFile(inFileName.Data());

  tree.SetEstimate(tree.GetEntries()+1);
  tree.Draw("num_LXeBb0n+num_LXeBb0n_eHi","covQual_sig == 3 && stat_sig == 0","goff");
  double median = TMath::Median(tree.GetSelectedRows(),tree.GetV1());
  double hl = getHalfLifeForCounts(yrs,median,xeMass);

  std::cout << "median counts: " << median << ", half-life = " << hl << std::endl;
  
  return hl;
}

//-----------------------------------------------------------
double nEXOUtil::getBkgdCounts(TString inFileName, int vol)
{
  // vol = 0: FV volume & full energy in SS
  // vol = 1: FV + FWHM
  // vol = 2: 3t + FWHM
  // vol = 3: 1t + FWHM
  
  TChain tree("tree");
  tree.AddFile(inFileName.Data());

  tree.SetEstimate(tree.GetEntries()+1);
  tree.Draw("bkg_tot:bkg_fwhm_fv:bkg_fwhm_3t:bkg_fwhm_1t:roi_bkg:roi_bkg_3t","covQual_sig == 3 && stat_sig == 0","para goff");

  double median = TMath::Median(tree.GetSelectedRows(),tree.GetVal(vol));

  //std::cout << "median bkg using vector " << vol << " = " << median << std::endl;

  return median;  
}


//-----------------------------------------------------------
Double_t nEXOUtil::getCutValue(TH1D* histo, Double_t gausSignif) {
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

//-----------------------------------------------------------
Double_t nEXOUtil::getHalfLifeForCounts(Double_t yrs, Double_t counts, Double_t xeMass) {
	return getNumXeAtoms(xeMass)*yrs*log(2.)/counts;
}

//-----------------------------------------------------------
Double_t nEXOUtil::getBkgdCountsInRange(const char* histFile, Double_t yrsExp, Int_t nPdfs, Double_t meanPerYear[], TString pdfNames[], Double_t* xRange, Double_t* yRange, Bool_t isBaTag) {
	TFile fHist(histFile);
	TH2D* h = 0;
	Double_t nBkgd = 0.;
	for (Int_t i=0; i<nPdfs; i++) {
		if (isBaTag && !pdfNames[i].Contains("bb2n")) continue;
		Double_t mean = meanPerYear[i];
		if (pdfNames[i].Contains("Co60")) {
			mean = getCo60CountsAfterYears(mean,yrsExp);
		} else {
			mean = mean*yrsExp;
		}
		h = (TH2D*)fHist.Get(Form("h_%s_ss", pdfNames[i].Data()));
		Double_t prob = getFractional2DIntegral(h, xRange[0], xRange[1], yRange[0], yRange[1]);

		nBkgd += prob*mean;
	}
	fHist.Close();
	return nBkgd;
}

//-----------------------------------------------------------
Double_t nEXOUtil::getCo60CountsAfterYears(Double_t meanPerYear, Double_t yrs) {
	// N0 = Ndecayed(1 year) / (1. - exp(-tau))
	// Ndecayed(years) = N0*(1.- exp(-tau*yrs))
	// tau = ln(2)/t_1/2, t_1/2 = 5.27yrs
	Double_t tau = log(2.)/5.27;
	Double_t scale = (1.-exp(-tau*yrs))/(1.-exp(-tau));
	return meanPerYear*scale;
}

//-----------------------------------------------------------
void nEXOUtil::getSmoothArray(Double_t* smoothArray, TH2* histo, Bool_t isSS) {
	Double_t t = (2457.83)/511.;
	Double_t me = 511.;
	Double_t p_ss[3] = {-2.244e-6, 36.90, 9.782e-3};
	Double_t p_ms[3] = {-3.676e-6, 41.08, 9.611e-3};

	//Get scale Factor for the resolution
	Double_t resAtQ = sqrt(p_ss[0]*p_ss[0]*t*me+p_ss[1]*p_ss[1]+p_ss[2]*p_ss[2]*t*me*t*me);
	Double_t resScaleFactor = 0.01*t*me/resAtQ;	

	//Make the energy observable
	RooRealVar energy("energy", "energy", 2100., 700., 3500.);
	energy.setBins(28000, "cache");
	RooArgSet obs(energy);

	//Make the Primakoff approx. pdf
	RooConstVar qValue("qValue", "qValue", t);
	TString bb2nFormula = "@0/511.*(@1-@0/511.)*(@1-@0/511.)*(@1-@0/511.)*(@1-@0/511.)*(@1-@0/511.)*(1.+2.*@0/511.+4./3.*@0/511.*@0/511.+@0/511.*@0/511.*@0/511./3.+@0/511.*@0/511.*@0/511.*@0/511./30.)";
	RooGenericPdf bb2nFunc("bb2nFunc", "Primakoff Approximation", bb2nFormula.Data(), RooArgList(energy, qValue));	
	
	//Make the Gaussian resolution pdf
	TString resFormula = Form("%03.f*sqrt(%0.3f*@0+%0.3f+%0.3f*@0*@0)", resScaleFactor, p_ss[0]*p_ss[0], p_ss[1]*p_ss[1], p_ss[2]*p_ss[2]);
	if (!isSS) 	resFormula = Form("%03.f*sqrt(%0.3f*@0+%0.3f+%0.3f*@0*@0)", resScaleFactor, p_ms[0]*p_ms[0], p_ms[1]*p_ms[1], p_ms[2]*p_ms[2]);
	RooFormulaVar width("width", "width", resFormula.Data(), RooArgList(energy));
	RooRealVar mean("mean", "mean", 0.);
	RooGaussian res("res", "res", energy, mean, width);

	//Make the convolution pdf, and extended convolution pdf
	RooFFTConvPdf	bb2nConv("bb2nConv", "bb2n", energy, bb2nFunc, res);
	bb2nConv.setBufferFraction(1.2);	
	RooRealVar norm("norm", "norm", 1.e7, 0., 1.e10);
	RooAddPdf bb2nExt("bb2nExt", "", RooArgSet(bb2nConv), RooArgSet(norm));

	//Get the dataset from the input histogram
	TH1D* hProjection = histo->ProjectionX("hProjection", 1, fNBinsY);
	RooDataHist data("data", "data", energy, hProjection);
	
	//Fit to the pdf in the specified range
	energy.setRange("r1", 2300., 3500.);
	if (!isSS) energy.setRange("r1", 2100., 3500.);
	bb2nExt.fitTo(data, RooFit::Range("r1"));

	//Get the integral in the fit range to calculate fitted normalization factor
	RooAbsReal* int1 = bb2nExt.createIntegral(energy, energy, "r1");
	Double_t i1 = int1->getVal();
	Double_t normFactor = norm.getVal()/i1/data.sumEntries();

	for (Int_t i=0; i<2800; i++) {
		energy.setVal(700.5+i*1.);
		smoothArray[i] = bb2nConv.getVal(&obs)*normFactor;
	}
}

//-----------------------------------------------------------
void nEXOUtil::smoothHisto(TH2* histo, Bool_t isSS, Int_t startBin) {
	//The following smooths out the bb2n pdf using the Primakoff approximation to the bb2n distribution
	Double_t smoothArray[fNBinsX];
	getSmoothArray(smoothArray, histo, isSS);

	Int_t nEntries = (Int_t)histo->Integral(1,fNBinsX,1,fNBinsY);

	Double_t totalVolume = TMath::Pi()*fFidVolXY*fFidVolXY*fFidVolZZ;

	Double_t s1[fNBinsY];
	Double_t s2[fNBinsY];
	Double_t shellVolume[fNBinsY];
	Double_t fraction[fNBinsY];
	Double_t shellCounts[fNBinsY];
	Double_t test = 0.;
	for (Int_t i=0; i<630; i++) {
		s1[i] = i*1.;
		s2[i] = (i+1)*1.;
		shellVolume[i] = TMath::Pi()*(fFidVolXY-s1[i])*(fFidVolXY-s1[i])*(fFidVolZZ-s1[i]);
		shellVolume[i] = shellVolume[i] - TMath::Pi()*(fFidVolXY-s2[i])*(fFidVolXY-s2[i])*(fFidVolZZ-s2[i]);
		fraction[i] = shellVolume[i]/totalVolume;
		shellCounts[i] = fraction[i]*nEntries;
		test+=fraction[i];
	}

	TH1D* hE = histo->ProjectionX("hE", 1, fNBinsY);

	for (Int_t i=startBin; i<=2000; i++) {
		hE->SetBinContent(i, smoothArray[i-1]*nEntries);
	}
	hE->Scale(1./hE->Integral());
		
	for (Int_t i=1; i<=630; i++) {
		for (Int_t j=1; j<=fNBinsX; j++) {
			Double_t binCon = hE->GetBinContent(j)*shellCounts[i-1];
			histo->SetBinContent(j, i, binCon);
		}
	}
	histo->Scale(1./histo->Integral(1,fNBinsX,1,fNBinsY)*nEntries);
}

//-----------------------------------------------------------
TCanvas* nEXOUtil::getEmptyDiscoveryPlot(enum MATRIX_ELEMENT matElem, Double_t xMin, Double_t xMax, Double_t yMin, Double_t yMax) {
	TString axisTitle = getMatElemAxisTitle(matElem);
	Double_t A = getMatElem(matElem);

	TString canvasName = Form("cc_%s", getMatElemName(matElem).Data());
	TString hFormatName = Form("hForm_%s", getMatElemName(matElem).Data());
	TGraph* gInv = getInvertedGraph(matElem);
	TGraph* gNor = getNormalGraph(matElem);

	TH1F* hFormat = new TH1F(hFormatName.Data(), "", 1, xMin, xMax);
	hFormat->GetXaxis()->CenterTitle();
	hFormat->GetXaxis()->SetTitleFont(22);
	hFormat->GetXaxis()->SetLabelFont(22);
	hFormat->GetXaxis()->SetTitle(fXaxisTitle.Data());
	hFormat->GetYaxis()->CenterTitle();
	hFormat->GetYaxis()->SetTitleFont(22);
	hFormat->GetYaxis()->SetLabelFont(22);
	hFormat->GetYaxis()->SetTitle(fYaxisTitle.Data());
	hFormat->GetYaxis()->SetTitleOffset(1.16);
	hFormat->GetYaxis()->SetRangeUser(yMin, yMax);
	hFormat->SetStats(false);

	TCanvas* cc = new TCanvas(canvasName, canvasName);
	hFormat->Draw();
	gInv->Draw("sameLF");
	gNor->Draw("sameLF");
	cc->Modified();
	cc->Update();		

	TGaxis *axis = new TGaxis(cc->GetUxmax(), cc->GetUymax(), cc->GetUxmax()*0.999, cc->GetUymin(), 
			1000.0*sqrt(A*1.0e24/yMax), 1000.0*sqrt(A*1.0e24/yMin), 510, "-LG");
	axis->SetTitleOffset(1.1);
	axis->SetLabelFont(22);
	axis->SetTitleFont(22);
	axis->CenterTitle();
	axis->SetTitle(axisTitle.Data());
	axis->Draw();

	cc->SetLogy();
	cc->SetTickx();
	cc->RedrawAxis();
	cc->Modified();
	cc->Update();

	return cc;
}

//-----------------------------------------------------------
TGraph* nEXOUtil::getInvertedGraph(enum MATRIX_ELEMENT matElem) {
	TString name = Form("g_inv_%s", getMatElemName(matElem).Data());
	Double_t A = getMatElem(matElem);

	Double_t invertedPoInt_ts_x[220];
	Double_t invertedPoInt_ts_y[220];
	for (Int_t i = 0; i < 110; i++){
		invertedPoInt_ts_x[i] = i*0.1;
		invertedPoInt_ts_y[i] = A*(1e24)/(0.015*0.015);
		invertedPoInt_ts_x[220-i-1] = i*0.1;
		invertedPoInt_ts_y[220-i-1] = A*(1e24)/(0.05*0.05);
	}

	TGraph* inverted = new TGraph(220, invertedPoInt_ts_x, invertedPoInt_ts_y);
	inverted->SetName(name.Data());
	inverted->SetTitle("");
	inverted->SetFillColor(kBlue-9);
	inverted->SetLineColor(kBlue-9);
	inverted->SetFillStyle(3002);

	return inverted;
}

//-----------------------------------------------------------
TGraph* nEXOUtil::getNormalGraph(enum MATRIX_ELEMENT matElem) {
	TString name = Form("g_inv_%s", getMatElemName(matElem).Data());
	Double_t A = getMatElem(matElem);

	//Now make the graph for the normal hierarchy
	Double_t normalPoInt_ts_x[220];
	Double_t normalPoInt_ts_y[220];
	for (Int_t i = 0; i < 110; i++){
		normalPoInt_ts_x[i] = i*0.1;
		normalPoInt_ts_y[i] = A*(1e24)/(0.0001*0.0001);
		normalPoInt_ts_x[220-i-1] = i*0.1;
		normalPoInt_ts_y[220-i-1] = A*(1e24)/(0.005*0.005);
	}
	TGraph* normal = new TGraph(220, normalPoInt_ts_x, normalPoInt_ts_y);
	normal->SetName(name.Data());
	normal->SetTitle("");
	normal->SetFillColor(kOrange-3);
	normal->SetLineColor(kOrange-3);
	normal->SetFillStyle(3002);

	return normal;
}

//-----------------------------------------------------------
std::vector<Double_t> nEXOUtil::getListOfDiscoveryCounts(std::vector<Double_t> bkgCounts, Double_t gausSignif, Double_t prob, Double_t eff) {
	std::vector<Double_t> retVec;
	Int_t size = (Int_t)bkgCounts.size();
	for (Int_t i=0; i<size; i++) {
		retVec.push_back(getDiscoveryCounts(bkgCounts.at(i), gausSignif, prob, eff));
	}
	return retVec;
}

//-----------------------------------------------------------
Double_t nEXOUtil::getDiscoveryCounts(Double_t bkgCounts, Double_t gausSignif, Double_t prob, Double_t eff) {
	Double_t alpha = ROOT::Math::normal_cdf_c(gausSignif, 1.);
	Double_t nCrit = getCriticalCounts(alpha, bkgCounts);
	Double_t nLds  = getLeastDetectableSignal(prob, bkgCounts, nCrit);
	return nLds/eff;
}

//-----------------------------------------------------------
Double_t nEXOUtil::getCriticalCounts(Double_t alpha, Double_t nBkgd) {
	//Takes in the required significance and background counts and returns the critical number of counts
	Double_t prob = 1.;
	Int_t counts = 0;
	while (prob > alpha) {
		prob = ROOT::Math::poisson_cdf_c(counts, nBkgd);
		counts++;
	}
	return (Double_t)counts;
}

//-----------------------------------------------------------
Double_t nEXOUtil::getLeastDetectableSignal(Double_t prob, Double_t nBkgd, Double_t nCrit) {
	//Takes in the required probability, background counts, and critical counts and returns the least detectable signal counts
	TF1 f("func", Form("ROOT::Math::poisson_cdf_c(%i, x)-%0.2f", (Int_t)nCrit, prob), 0., 1000.);
	ROOT::Math::WrappedTF1 wf1(f);

	ROOT::Math::BrentRootFinder brf;
	brf.SetNpx(50000);
	brf.SetFunction(wf1, 0., 1000.);
	brf.Solve();

	return brf.Root()-nBkgd;
}

//-----------------------------------------------------------
void nEXOUtil::normalizeHistoFile(const char* inFileName, const char* outFileName) {
	//Takes in a name for a ROOT file containing TH1 histograms, creates a normalized copy and saves to outFileName
	TFile fIn(inFileName);
	TList* list = fIn.GetListOfKeys();

	TFile fOut(outFileName, "RECREATE");
	TH1* h = 0;
	for (Int_t nn=0; nn<list->GetSize(); nn++) {
		TString name = list->At(nn)->GetName();
		h = (TH1*)fIn.Get(name.Data());
		h->SetName(Form("%s_old", name.Data()));
		TH1* hnew = (TH1*)h->Clone(name.Data());
		hnew->Scale(1./hnew->Integral());
		fOut.cd();
		hnew->Write();
	}
	fOut.Close();
	fIn.Close();
}

//-----------------------------------------------------------
void nEXOUtil::rebin2DHistoFile(const char* inFileName, const char* outFileName, Int_t rebinX, Int_t rebinY) {
	//Rebins the histograms in inFileName given the number of bins to combine Int_to 1 (all bin widths are assumed to be 1 keV or 1mm to start with). It also removes standoff bins above 640mm since these are empty by definition of the standoff-distance.
	TFile fIn(inFileName);
	TList* list = fIn.GetListOfKeys();

	TFile fOut(outFileName, "RECREATE");
	TH2D* h = 0;
	for (Int_t nn=0; nn<list->GetSize(); nn++) {
		TString name = list->At(nn)->GetName();
		h = (TH2D*)fIn.Get(name.Data());
		h->SetName(Form("%s_old", name.Data()));
		TH2D* hnew = new TH2D(name.Data(), "", 2800, 700., 3500., 640, 0., 640.);
		for (Int_t i=1; i<=2800; i++) {
			for (Int_t j=1; j<=640; j++) {
					hnew->SetBinContent(i,j,h->GetBinContent(i,j));
			}
		}
		hnew->RebinX(rebinX);
		hnew->RebinY(rebinY);
		fOut.cd();
		hnew->Write();
	}
	fOut.Close();
	fIn.Close();
}

//-----------------------------------------------------------
TH2D* nEXOUtil::combineHistos(std::vector<TH2D*>* histos, std::vector<Double_t>* fracs, Double_t norm) {
	//this script takes in a vector of TH1 objects and normalization fractions and returns a combined sum TH1 according to the given fractions. Optional overall normalization can be set.

	//Do some consistency checks
	if (histos->size() != fracs->size()) {
		std::cout<<"Input vectors must contain the same number of entries!"<<std::endl;
		return NULL;
	}
	Int_t size = (Int_t)histos->size();

	Int_t nBinsX = histos->at(0)->GetXaxis()->GetNbins();
	Int_t nBinsY = histos->at(0)->GetYaxis()->GetNbins();
	Double_t xlo = histos->at(0)->GetXaxis()->GetXmin();
	Double_t xhi = histos->at(0)->GetXaxis()->GetXmax();
	Double_t ylo = histos->at(0)->GetYaxis()->GetXmin();
	Double_t yhi = histos->at(0)->GetYaxis()->GetXmax();

	TH2D* retHist = new TH2D("retHist", "", nBinsX, xlo, xhi, nBinsY, ylo, yhi);
	std::cout << retHist->GetNbinsX() << std::endl;
	std::cout << retHist->GetNbinsY() << std::endl;
	for (Int_t i=0; i<size; i++) {
		std::cout << histos->at(i)->GetName() << " , " << histos->at(i)->GetNbinsX() << " , " << histos->at(i)->GetNbinsY() << std::endl;	
		retHist->Add(histos->at(i), fracs->at(i));
	}
	retHist->Scale(norm/retHist->Integral());
	return retHist;
}

//-----------------------------------------------------------
Double_t nEXOUtil::getFractional2DIntegral(TH2* h, Double_t xLo, Double_t xHi, Double_t yLo, Double_t yHi) {
	//this script takes in a 2D histogram and returns the Int_tegral in the given range with constant Int_terpolation used if the range extrema don't lie on bin edges between the bins

	Int_t    xBins  = h->GetXaxis()->GetNbins(), 		yBins  = h->GetYaxis()->GetNbins();
	Double_t xWidth = h->GetXaxis()->GetBinWidth(1), 	yWidth = h->GetYaxis()->GetBinWidth(1);

	Int_t xLoBin = h->GetXaxis()->FindBin(xLo), yLoBin = h->GetYaxis()->FindBin(yLo);
	Int_t xHiBin = h->GetXaxis()->FindBin(xHi), yHiBin = h->GetYaxis()->FindBin(yHi);
	if (xHiBin>xBins) xHiBin=xBins;
	if (yHiBin>yBins) yHiBin=yBins;
	if (xLoBin<=0) xLoBin=1;
	if (yLoBin<=0) yLoBin=1;
	
	//Get the fractions of the bins contained in the Int_tegral range
	Double_t lFrac = ( h->GetXaxis()->GetBinUpEdge (xLoBin)-xLo )/xWidth;
	Double_t rFrac = (-h->GetXaxis()->GetBinLowEdge(xHiBin)+xHi )/xWidth;
	Double_t bFrac = ( h->GetYaxis()->GetBinUpEdge (yLoBin)-yLo )/yWidth;
	Double_t tFrac = (-h->GetYaxis()->GetBinLowEdge(yHiBin)+yHi )/yWidth;

	//Calculate the Int_tegral inside the range (not includeing fractional bins)
	Double_t inInt = h->Integral(xLoBin+1, xHiBin-1, yLoBin+1, yHiBin-1);

	//Calculate the Int_tegral of the edges corresponding to the Int_tegral range, excluding corners
	Double_t lInt = h->Integral(xLoBin, xLoBin, yLoBin+1, yHiBin-1)*lFrac;
	Double_t rInt = h->Integral(xHiBin, xHiBin, yLoBin+1, yHiBin-1)*rFrac;
	Double_t bInt = h->Integral(xLoBin+1, xHiBin-1, yLoBin, yLoBin)*bFrac;
	Double_t tInt = h->Integral(xLoBin+1, xHiBin-1, yHiBin, yHiBin)*tFrac;

	//Calculate the Int_tegral of the four corners of the Int_tegral range
	Double_t lbInt = h->Integral(xLoBin, xLoBin, yLoBin, yLoBin)*lFrac*bFrac;
	Double_t ltInt = h->Integral(xLoBin, xLoBin, yHiBin, yHiBin)*lFrac*tFrac;
	Double_t rbInt = h->Integral(xHiBin, xHiBin, yLoBin, yLoBin)*rFrac*bFrac;
	Double_t rtInt = h->Integral(xHiBin, xHiBin, yHiBin, yHiBin)*rFrac*tFrac;
		
	//return the Int_tegral
	return inInt + lInt + rInt + bInt + tInt + lbInt + ltInt + rbInt + rtInt;
}
