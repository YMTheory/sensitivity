#include "nEXOMbbVsMnuPlot.hh"

ClassImp(nEXOMbbVsMnuPlot);

nEXOMbbVsMnuPlot::nEXOMbbVsMnuPlot(const char* name, const char* title)
  : nEXOSensPlot(name,title), fNuMassOption(kNuMassMin), fNpoints(1000),
    fLoXRangeCoeff(1.), fLoXRangeExp(3), fHiXRangeCoeff(1.), fHiXRangeExp(0),
    fLoYRangeCoeff(1.), fLoYRangeExp(3), fHiYRangeCoeff(1.), fHiYRangeExp(0)
{
  CreateAvailableBands();
  nEXONuOscPars::GetInstance()->SetFitSource(nEXONuOscPars::kprd90_093006_2014);
}

nEXOMbbVsMnuPlot::~nEXOMbbVsMnuPlot() {}

TObject* nEXOMbbVsMnuPlot::GetPlot()
{
  // main function to return plot of half-life vs detector livetimembb vs nu mass (min, sum, etc)

  SetBand("exo200_current","EXO200, Current",0.190,0.450);
  SetBand("exo_ultimate","EXO200, Ultimate",0.05,0.1334);

  //SetBand("nexo_sens_gamma","nEXO, No Ba-tagging",0.0072, 0.0192);
  //SetBand("nexo_sens_batag","nEXO, With Ba-tagging",0.0031, 0.0083);

  
  TCanvas* canvas = CreateEmptyCanvas();
  
  for(std::map<TString, TGraph>::iterator band = fBands.begin(); band != fBands.end(); band++)
  {
    TString bandName = band->first;
    std::cout << "Plotting band " << bandName << std::endl;
    canvas->cd();
    TGraph* g_band = &band->second;
    g_band->Draw("same f");
    g_band->Draw("same l");
    TPaveText* p_band = &fTitles.at(bandName);
    p_band->Draw();    
  }
  

  gStyle->SetLineWidth(3);
  gPad->RedrawAxis();

  canvas->cd();
  canvas->Modified();
  canvas->Update();

  
  return canvas;  
}

bool nEXOMbbVsMnuPlot::SetBand(const char* name, const char* title, double min, double max)
{
  TString bandName(name);
  if(fLineColors.count(bandName.Data()) <= 0)
    return false;
  TString bandTitle(title);

  std::cout << "Setting band " << bandName << " : " << bandTitle << " with " << min << " , " << max << std::endl;
    
  Double_t loRange =  0.99*GetLoXRange();
  Double_t hiRange =  1.01*GetHiXRange();
  Double_t limits[4] = {hiRange, loRange, loRange, hiRange};

  Double_t points[4] = {min, min, max, max};
  TGraph graph(4, limits, points);

  fBands.insert(std::make_pair(bandName,TGraph(graph)));
  fBands.at(bandName).SetName(bandName.Data());
  fBands.at(bandName).SetTitle(bandTitle.Data());

  TGraph *band = &fBands.at(bandName);
  
  band->SetLineColor(fLineColors.at(bandName));
  band->SetLineWidth(fLineWidths.at(bandName));
  band->SetFillColor(fFillColors.at(bandName));
  band->SetFillStyle(fFillStyles.at(bandName));  

  Double_t xCenter = sqrt(GetHiXRange() * GetLoXRange());
  Double_t yCenter = sqrt(min * max);

  TPaveText pave(xCenter,yCenter,xCenter,yCenter, "br");
  pave.SetBorderSize(0);
  pave.SetFillColor(0);
  pave.SetFillStyle(0);
  TText* paveText = pave.AddText(bandTitle.Data());
  paveText->SetTextColor(fTextColors.at(bandName));
  paveText->SetTextFont(fTextFonts.at(bandName));
  paveText->SetTextSize(fTextSizes.at(bandName));
  
  fTitles.insert(std::make_pair(bandName,TPaveText(pave)));
     
  return true;
}

bool nEXOMbbVsMnuPlot::SetBandFromFile(const char* name, const char* title, const char* filename, double yrs, double mass)
{
  TString bandName(name);
  if(fLineColors.count(bandName.Data()) <= 0)
    return false;
  TString bandTitle(title);

  Double_t hl = nEXOUtils::GetSensHalfLife(filename,yrs,mass);
  double min = 1e6;
  double max = -1e6;
  for(size_t i = 0; i < nEXONuclearMatrixElement::nNMEs; i++)
  {
    nEXONuclearMatrixElement nme(static_cast<nEXONuclearMatrixElement::NME_t>(i));
    double m = sqrt(nme.GetAxisScale()*1e24/hl);
    if(m < min)
      min = m;
    if(m > max)
      max = m;
  }

  return SetBand(bandName.Data(),bandTitle.Data(),min,max);  
}

TCanvas* nEXOMbbVsMnuPlot::CreateEmptyCanvas()
{
  TH2D* h_axis = new TH2D(Form("h_axis_%s_%d",GetName(),fNuMassOption), "", fNpoints, GetLoXRange(), GetHiXRange(), fNpoints, GetLoYRange(), GetHiYRange());
  h_axis->SetStats(false);

  TString xTitle = "m_{min} (eV)";
  if(fNuMassOption == kNuMassSum)
    xTitle = "#Sigmam_{i} (eV)";
  h_axis->GetXaxis()->SetTitle(xTitle.Data());
  
  h_axis->GetYaxis()->SetTitle("m_{#beta#beta} (eV)");
  h_axis->GetXaxis()->CenterTitle();
  h_axis->GetXaxis()->SetTitleSize(0.04);	
  h_axis->GetYaxis()->CenterTitle();
  h_axis->GetYaxis()->SetTitleSize(0.04);
  h_axis->GetXaxis()->SetLabelSize(0.04);
  h_axis->GetYaxis()->SetLabelSize(0.04);
  
  TGraph* g_i_cv = new TGraph(fNpoints*2);
  TGraph* g_n_cv = new TGraph(fNpoints*2);
  TGraph* g_i_err = new TGraph(fNpoints*2);
  TGraph* g_n_err = new TGraph(fNpoints*2);

  g_i_cv->SetFillColor(kRed-9);
  g_n_cv->SetFillColor(kBlue-5);
  g_i_err->SetFillColor(0);
  g_n_err->SetFillColor(0);
  
  g_i_err->SetLineWidth(2);
  g_n_err->SetLineWidth(2);

  g_i_cv->SetFillStyle(1001);
  g_n_cv->SetFillStyle(1001);
  g_i_err->SetFillStyle(0);
  g_n_err->SetFillStyle(0);
  
  g_i_cv->SetLineColor(kRed);
  g_n_cv->SetLineColor(kBlue);
  g_i_err->SetLineColor(kRed);
  g_n_err->SetLineColor(kBlue);
  
  g_i_err->SetLineStyle(kDashed);
  g_n_err->SetLineStyle(kDashed);

  switch(fNuMassOption)
  {
      case kNuMassMin:
        FillGraphPointsMmin(g_n_cv,g_n_err,g_i_cv,g_i_err);
        break;
      case kNuMassSum:
        FillGraphPointsMsum(g_n_cv,g_n_err,g_i_cv,g_i_err);
  }
  
  TCanvas* cc = new TCanvas(Form("%s_%d",GetName(),fNuMassOption),GetTitle());//"", 3000, 3000);
  cc->SetTickx();
  cc->SetTicky();
  cc->SetLogy();
  cc->SetLogx();
  h_axis->Draw();
  
  g_n_cv->Draw("same f");
  g_n_cv->Draw("same l");
  g_n_err->Draw("same f");
  g_n_err->Draw("same l");
  g_i_cv->Draw("same f");
  g_i_cv->Draw("same l");
  g_i_err->Draw("same f");
  g_i_err->Draw("same l");
  
  cc->RedrawAxis();
  return cc;
}

void nEXOMbbVsMnuPlot::FillGraphPointsMmin(TGraph* g_n_cv,TGraph* g_n_err,TGraph* g_i_cv,TGraph* g_i_err)
{
  Int_t nLo = 0;//, nLo_68 = 0, nLo_90 = 0, nLo_95 = 0;
  Double_t mbb[4];
  Double_t mbb_lo[4];
  Double_t mbb_hi[4];
  Double_t first = TMath::Log10(fLoXRangeCoeff) - fLoXRangeExp;
  Double_t last = TMath::Log10(fHiXRangeCoeff) + fHiXRangeExp;
  Double_t step = (last - first)*1./(fNpoints-1);
  
  for (Int_t i=0; i<fNpoints; i++) {
    //Double_t exponent = fLoXRangeExp*(-1. + (Double_t)i/(fNpoints-1));
    //Double_t m_min = fLoXRangeCoeff*TMath::Power(10., exponent);
    Double_t exponent = first + i*step;
    Double_t m_min = TMath::Power(10., exponent);
    
    //Make the graphs for the central values
    nEXONuOscPars::GetInstance()->EvalNormalMbbMmin(mbb, mbb_lo, mbb_hi, m_min, 2); // eval all possible mbbs from nu osc limits

    // find actual min and max mbb
    Int_t minBin = 0, maxBin = 0;
    Double_t min = 1.e6, min_lo=0.;
    Double_t max = -1.e6, max_hi=0.;
    for (int j=0; j<4; j++) {
      if (mbb[j] < min) {min = mbb[j]; minBin = j;}
      if (mbb[j] > max) {max = mbb[j]; maxBin = j;}
    }
    // switch with axis limit if outside of it
    if (min < GetLoYRange() && nLo==0) nLo = 1;
    if (min > GetLoYRange() && nLo==1) nLo = 2;
    if (min < GetLoYRange() && nLo==2) nLo = 3;
    if (min > GetLoYRange() && nLo==3) nLo = 4;		
    if (nLo==1 || nLo==2 || nLo==3) min = GetLoYRange();
    min_lo = min-mbb_lo[minBin];
    max_hi = max+mbb_hi[maxBin];
    
    g_n_cv->SetPoint(i, m_min, max);
    g_n_cv->SetPoint(2*fNpoints-i-1, m_min, min);
    g_n_err->SetPoint(i, m_min, max_hi);
    g_n_err->SetPoint(2*fNpoints-i-1, m_min, min_lo);	
    
    nEXONuOscPars::GetInstance()->EvalInvertedMbbMmin(mbb, mbb_lo, mbb_hi, m_min, 2);
    minBin = 0, maxBin = 0;
    min = 1.e6, min_lo=0.;
    max = -1.e6, max_hi=0.;
    for (int j=0; j<4; j++) {
      if (mbb[j] < min) {min = mbb[j]; minBin = j;}
      if (mbb[j] > max) {max = mbb[j]; maxBin = j;}
    }
    min_lo = min-mbb_lo[minBin];
    max_hi = max+mbb_hi[maxBin];
    
    g_i_cv->SetPoint(i, m_min, max);
    g_i_cv->SetPoint(2*fNpoints-i-1, m_min, min);	
    g_i_err->SetPoint(i, m_min, max_hi);
    g_i_err->SetPoint(2*fNpoints-i-1, m_min, min_lo);
  }
}


void nEXOMbbVsMnuPlot::FillGraphPointsMsum(TGraph* g_n_cv,TGraph* g_n_err,TGraph* g_i_cv,TGraph* g_i_err)
{
  Int_t nLo = 0;//, nLo_68 = 0, nLo_90 = 0, nLo_95 = 0;
  Double_t mbb[4];
  Double_t mbb_lo[4];
  Double_t mbb_hi[4];
  Double_t first = TMath::Log10(fLoXRangeCoeff) - fLoXRangeExp;
  Double_t last = TMath::Log10(fHiXRangeCoeff) + fHiXRangeExp;
  Double_t step = (last - first)*1./(fNpoints-1);

  Int_t nBlankPointsNormal = 0;
  Int_t nBlankPointsInverted = 0;
  
  for (Int_t i=0; i<fNpoints; i++) {
    Double_t exponent = first + i*step;
    Double_t m_sum = TMath::Power(10., exponent);//Double_t m_min = TMath::Power(10., exponent);
    
    //Make the graphs for the central values
    if(not nEXONuOscPars::GetInstance()->EvalNormalMbbMsum(mbb, mbb_lo, mbb_hi, m_sum, 2)) // eval all possible mbbs from nu osc limits
    {
      nBlankPointsNormal++;
      continue;
    }
      
    // find actual min and max mbb
    Int_t minBin = 0, maxBin = 0;
    Double_t min = 1.e6, min_lo=0.;
    Double_t max = -1.e6, max_hi=0.;
    for (int j=0; j<4; j++) {
      if (mbb[j] < min) {min = mbb[j]; minBin = j;}
      if (mbb[j] > max) {max = mbb[j]; maxBin = j;}
    }
    // switch with axis limit if outside of it
    if (min < GetLoYRange() && nLo==0) nLo = 1;
    if (min > GetLoYRange() && nLo==1) nLo = 2;
    if (min < GetLoYRange() && nLo==2) nLo = 3;
    if (min > GetLoYRange() && nLo==3) nLo = 4;		
    if (nLo==1 || nLo==2 || nLo==3) min = GetLoYRange();
    min_lo = min-mbb_lo[minBin];
    max_hi = max+mbb_hi[maxBin];
    
    g_n_cv->SetPoint(i, m_sum, max);
    g_n_cv->SetPoint(2*fNpoints-i-1, m_sum, min);
    g_n_err->SetPoint(i, m_sum, max_hi);
    g_n_err->SetPoint(2*fNpoints-i-1, m_sum, min_lo);
  }

  for (Int_t i=0; i<fNpoints; i++) {
    Double_t exponent = first + i*step;
    Double_t m_sum = TMath::Power(10., exponent);//Double_t m_min = TMath::Power(10., exponent);

    if(not nEXONuOscPars::GetInstance()->EvalInvertedMbbMsum(mbb, mbb_lo, mbb_hi, m_sum, 2)) // eval all possible mbbs from nu osc limits
    {
      nBlankPointsInverted++;
      continue;
    }
      
    // find actual min and max mbb
    Int_t minBin = 0, maxBin = 0;
    Double_t min = 1.e6, min_lo=0.;
    Double_t max = -1.e6, max_hi=0.;
    for (int j=0; j<4; j++) {
      if (mbb[j] < min) {min = mbb[j]; minBin = j;}
      if (mbb[j] > max) {max = mbb[j]; maxBin = j;}
    }
    min_lo = min-mbb_lo[minBin];
    max_hi = max+mbb_hi[maxBin];
    
    g_i_cv->SetPoint(i, m_sum, max);
    g_i_cv->SetPoint(2*fNpoints-i-1, m_sum, min);	
    g_i_err->SetPoint(i, m_sum, max_hi);
    g_i_err->SetPoint(2*fNpoints-i-1, m_sum, min_lo);
    
  }

  for(Int_t i=0; i<nBlankPointsNormal; i++)
  {
    g_n_cv->RemovePoint(0);
    g_n_cv->RemovePoint(g_n_cv->GetN()-1);
    g_n_err->RemovePoint(0);
    g_n_err->RemovePoint(g_n_err->GetN()-1);
  }
  for(Int_t i=0; i<nBlankPointsInverted; i++)
  {
    g_i_cv->RemovePoint(0);
    g_i_cv->RemovePoint(g_i_cv->GetN()-1);
    g_i_err->RemovePoint(0);
    g_i_err->RemovePoint(g_i_err->GetN()-1);
  }
}


void nEXOMbbVsMnuPlot::PlotEXO200(TCanvas& canvas)
{
  Double_t loRange =  0.99*GetLoXRange();
  Double_t hiRange =  1.01*GetHiXRange();
  Double_t limits[4] = {hiRange, loRange, loRange, hiRange};

  //Double_t exo200Current[4] 	= {0.1392, 0.1392, 0.3717, 0.3717}; //2012 result
  Double_t exo200CurrentLo = 0.190;
  Double_t exo200CurrentHi = 0.450;
  Double_t exo200Current[4]  = {exo200CurrentLo, exo200CurrentLo, exo200CurrentHi, exo200CurrentHi};

  Double_t exo200UltimateLo = 0.0500;
  Double_t exo200UltimateHi = 0.1334;    
  Double_t exo200Ultimate[4] = {exo200UltimateLo, exo200UltimateLo, exo200UltimateHi, exo200UltimateHi};
  TGraph* g_exo200Current 	= new TGraph(4, limits, exo200Current);
  TGraph* g_exo200Ultimate 	= new TGraph(4, limits, exo200Ultimate);
  
  g_exo200Current->SetLineWidth(2);
  g_exo200Current->SetLineColor(kBlue+1);
  g_exo200Current->SetFillColor(kBlue-9);
  g_exo200Current->SetFillStyle(3001);
  
  g_exo200Ultimate->SetLineWidth(2);
  g_exo200Ultimate->SetLineColor(kMagenta+1);
  g_exo200Ultimate->SetFillColor(kMagenta-9);
  g_exo200Ultimate->SetFillStyle(3001);
  
  canvas.cd();
  g_exo200Current->Draw("same f");
  g_exo200Current->Draw("same l");
  g_exo200Ultimate->Draw("same f");
  g_exo200Ultimate->Draw("same l");

  // center is given by geometric average pow^1/2 as oppose to sum-1/2 because of the axis log-scale
  Double_t xCenter = sqrt(GetHiXRange() * GetLoXRange());
  Double_t yCenter = sqrt(exo200CurrentLo * exo200CurrentHi);
  TPaveText* exo200Current_title = new TPaveText(xCenter,yCenter,xCenter,yCenter, "br"); //TPaveText(0.5, 0.76, 0.5, 0.76, "brNDC");
  exo200Current_title->SetBorderSize(0);
  exo200Current_title->SetFillColor(0);
  exo200Current_title->SetFillStyle(0);
  TText* exo200Current_title_text = exo200Current_title->AddText("EXO200, Current");
  exo200Current_title_text->SetTextColor(kBlack);
  exo200Current_title_text->SetTextFont(22);
  exo200Current_title_text->SetTextSize(0.05);

  yCenter = sqrt(exo200UltimateLo * exo200UltimateHi);
  TPaveText* exo200Ultimate_title = new TPaveText(xCenter,yCenter,xCenter,yCenter, "br");
  exo200Ultimate_title->SetBorderSize(0);
  exo200Ultimate_title->SetFillColor(0);
  exo200Ultimate_title->SetFillStyle(0);
  TText* exo200Ultimate_title_text = exo200Ultimate_title->AddText("EXO200, Ultimate");
  exo200Ultimate_title_text->SetTextColor(kBlack);
  exo200Ultimate_title_text->SetTextFont(22);
  exo200Ultimate_title_text->SetTextSize(0.05);
  
  exo200Current_title->Draw();
  exo200Ultimate_title->Draw();

  return ;
}

void nEXOMbbVsMnuPlot::CreateAvailableBands()
{

  TString name;

  name = "exo200_current";  
  fLineColors[name] = kBlue+1; fLineWidths[name] = 2; fFillColors[name] = kBlue-9;  fFillStyles[name] = 3001;
  fTextColors[name] = kBlack; fTextFonts[name] = 22; fTextSizes[name] = 0.05;

  name = "exo_ultimate";
  fLineColors[name] = kMagenta+1; fLineWidths[name] = 2; fFillColors[name] = kMagenta-9;  fFillStyles[name] = 3001;
  fTextColors[name] = kBlack; fTextFonts[name] = 22; fTextSizes[name] = 0.05;

  name = "nexo_sens_gamma";
  fLineColors[name] = kGreen+1; fLineWidths[name] = 2; fFillColors[name] = kGreen-9;  fFillStyles[name] = 3001;
  fTextColors[name] = kBlack; fTextFonts[name] = 22; fTextSizes[name] = 0.05;

  name = "nexo_sens_batag";
  fLineColors[name] = kRed+1; fLineWidths[name] = 2; fFillColors[name] = kRed-9;  fFillStyles[name] = 3001;
  fTextColors[name] = kBlack; fTextFonts[name] = 22; fTextSizes[name] = 0.05;
  
}
