#include "nEXOHalflifeVsTimePlot.hh"

ClassImp(nEXOHalflifeVsTimePlot);

nEXOHalflifeVsTimePlot::nEXOHalflifeVsTimePlot(const char* name, const char* title) : nEXOSensPlot(name,title), fNME(0)
{
  fXaxisTitle = "Exposure (yrs)";
  fYaxisTitle = "Xe-136 #beta#beta0#nu T_{1/2} (yrs)";//"^{136}Xe T_{1/2}";
  
  SetAxisLimits();

  CreateAvailableGraphs();
}

nEXOHalflifeVsTimePlot::~nEXOHalflifeVsTimePlot() {}

void nEXOHalflifeVsTimePlot::SetAxisLimits(Double_t xMin, Double_t xMax, Double_t yMin, Double_t yMax)
{
  fXmin = xMin; fXmax = xMax;
  fYmin = yMin; fYmax = yMax;
}

void nEXOHalflifeVsTimePlot::CreateAvailableGraphs()
{
  fLineColors["sens_90"] = kBlue;
  fLineWidths["sens_90"] = 2;
  fLineStyles["sens_90"] = 1;
  fSmoothFunction["sens_90"] = new TF1("f_sens_90","[0]*pow(x,[1])");
  fSmoothFunction["sens_90"]->SetParameters(1,1);

  fLineColors["disc3_50"] = kRed;
  fLineWidths["disc3_50"] = 2;
  fLineStyles["disc3_50"] = 1;
  fSmoothFunction["disc3_50"] = new TF1("f_disc3_50","[0]*pow(x,[1])");
  fSmoothFunction["disc3_50"]->SetParameters(1,1);

  fLineColors["count_sens_90"] = kBlue-9;
  fLineWidths["count_sens_90"] = 2;
  fLineStyles["count_sens_90"] = 2;
  fSmoothFunction["count_sens_90"] = new TF1("f_count_sens_90","[0]*pow(x,[1])");
  fSmoothFunction["count_sens_90"]->SetParameters(1,1);

  fLineColors["count_disc3_50"] = kRed-9;
  fLineWidths["count_disc3_50"] = 2;
  fLineStyles["count_disc3_50"] = 2;
  fSmoothFunction["count_disc3_50"] = new TF1("f_disc3_50","[0]*pow(x,[1])");
  fSmoothFunction["count_disc3_50"]->SetParameters(1,1);
}


void nEXOHalflifeVsTimePlot::SetNME(nEXONuclearMatrixElement::NME_t nme)
{
  if(fNME)
    delete fNME;
  fNME = new nEXONuclearMatrixElement(nme);
}

TObject* nEXOHalflifeVsTimePlot::GetPlot()
{
  // main function to return plot of half-life vs detector livetime

  TCanvas* canvas = CreateEmptyCanvas();
  TLegend* leg = new TLegend(0.47, 0.13, 0.87, 0.28);
  leg->SetMargin(0.13);

  for(std::map<TString, TGraph>::iterator graph = fGraphs.begin(); graph != fGraphs.end(); graph++)
    PlotGraph(graph->second,*canvas,*leg,fSmoothFunction.at(graph->first));
  
  PlotEXO200(*canvas,*leg);

  canvas->cd();
  leg->Draw();
  canvas->Modified();
  canvas->Update();
  
  return canvas;
}


bool nEXOHalflifeVsTimePlot::SetGraphPoints(const char* name, const char* title, size_t n, double* yrs, double* hls, double unit)
{
  TString graphName(name);

  if(fLineColors.count(graphName.Data()) <= 0)
    return false;

  for(size_t i = 0; i < n; i++)
    hls[i] *= unit;
  
  TGraph graph(n,yrs,hls);
  AddGraph(graphName, title, graph);
  
  return true;
}

void nEXOHalflifeVsTimePlot::AddGraph(TString name, TString title, TGraph& graph)
{
  fGraphs.insert(std::make_pair(name,TGraph(graph)));
  fGraphs.at(name).SetTitle(title.Data());
  SetGraphProperties(name.Data());
}

void nEXOHalflifeVsTimePlot::SetGraphProperties(const char* name)
{
  TGraph *graph = &fGraphs.at(name);
  
  graph->SetLineColor(fLineColors.at(name));
  graph->SetLineStyle(fLineStyles.at(name));
  graph->SetLineWidth(fLineWidths.at(name));
}


TGraph* nEXOHalflifeVsTimePlot::GetInvertedGraph(nEXONuclearMatrixElement& nme)
{
  // return graph with inverted neutrino mass ordering based with halflife conversion to mbb using nme
  
  TString name = Form("g_inv_%s", nme.GetName());
  Double_t A = nme.GetAxisScale();

  int max = static_cast<int>(11*fXmax);
  Double_t invertedPoInt_ts_x[2*max];
  Double_t invertedPoInt_ts_y[2*max];
  for (Int_t i = 0; i < max; i++){
    invertedPoInt_ts_x[i] = i*0.1;
    invertedPoInt_ts_y[i] = A*(1e24)/(0.015*0.015);
    invertedPoInt_ts_x[2*max-i-1] = i*0.1;
    invertedPoInt_ts_y[2*max-i-1] = A*(1e24)/(0.05*0.05);
  }

  TGraph* inverted = new TGraph(2*max, invertedPoInt_ts_x, invertedPoInt_ts_y);
  inverted->SetName(name.Data());
  inverted->SetTitle("");
  inverted->SetFillColor(kBlue-9);
  inverted->SetLineColor(kBlue-9);
  inverted->SetFillStyle(3002);
  
  return inverted;
}

TGraph* nEXOHalflifeVsTimePlot::GetNormalGraph(nEXONuclearMatrixElement& nme)
{
  // return graph with normal neutrino mass ordering based with halflife conversion to mbb using nme
  
  TString name = Form("g_inv_%s", nme.GetName());
  Double_t A = nme.GetAxisScale();

  int max = static_cast<int>(11*fXmax);
  Double_t normalPoInt_ts_x[2*max];
  Double_t normalPoInt_ts_y[2*max];
  for (Int_t i = 0; i < max; i++){
    normalPoInt_ts_x[i] = i*0.1;
    normalPoInt_ts_y[i] = A*(1e24)/(0.0001*0.0001);
    normalPoInt_ts_x[2*max-i-1] = i*0.1;
    normalPoInt_ts_y[2*max-i-1] = A*(1e24)/(0.005*0.005);
  }
  TGraph* normal = new TGraph(2*max, normalPoInt_ts_x, normalPoInt_ts_y);
  normal->SetName(name.Data());
  normal->SetTitle("");
  normal->SetFillColor(kOrange-3);
  normal->SetLineColor(kOrange-3);
  normal->SetFillStyle(3002);
  
  return normal;
}

TCanvas* nEXOHalflifeVsTimePlot::CreateEmptyCanvas()
{
  TString canvasName(GetName());

  if(fNME)
    canvasName = Form("%s_%s",canvasName.Data(),fNME->GetName());

  TCanvas* cc = new TCanvas(canvasName, GetTitle());

  TString hFormatName(Form("hForm_%s",canvasName.Data()));
  
  TH1D* hFormat = new TH1D(hFormatName.Data(),"", 1, fXmin, fXmax);
  hFormat->GetXaxis()->CenterTitle();
  hFormat->GetXaxis()->SetTitleFont(22);
  hFormat->GetXaxis()->SetTitleSize(0.06);
  hFormat->GetXaxis()->SetTitleOffset(0.75);
  hFormat->GetXaxis()->SetLabelFont(22);
  hFormat->GetXaxis()->SetLabelSize(0.05);
  hFormat->GetXaxis()->SetLabelOffset(0.);
  hFormat->GetXaxis()->SetTitle(fXaxisTitle.Data());
  hFormat->GetYaxis()->CenterTitle();
  hFormat->GetYaxis()->SetTitleFont(22);
  hFormat->GetYaxis()->SetTitleSize(0.06);
  hFormat->GetYaxis()->SetTitleOffset(0.8);
  hFormat->GetYaxis()->SetLabelFont(22);
  hFormat->GetYaxis()->SetLabelSize(0.05);
  hFormat->GetYaxis()->SetLabelOffset(0.);
  hFormat->GetYaxis()->SetTitle(fYaxisTitle.Data());
  hFormat->GetYaxis()->SetRangeUser(fYmin,fYmax);
  hFormat->SetStats(false);

  cc->cd();
  hFormat->Draw();

  if(fNME)
  {
    Double_t A = fNME->GetAxisScale();
    TString axisTitle = fNME->GetTitle();

    hFormat->SetTitle(axisTitle.Data());

    TGraph* gInv = GetInvertedGraph(*fNME);
    TGraph* gNor = GetNormalGraph(*fNME);
    
    gInv->Draw("sameF");
    gNor->Draw("sameF");

    Double_t xCenter = (fXmax + fXmin)/2.;
    //Double_t yCenter = FindMinGraph(gNor);//sqrt(fYmin * fYmax);
    Double_t xmin, ymin, xmax, ymax;
    gInv->ComputeRange(xmin,ymin,xmax,ymax);
    Double_t yCenter = sqrt(ymin * ymax);//A*(1e24)/(0.005*0.005);
    //std::cout << axisTitle << " A " << A << " = " << yCenter << std::endl;
 
    TPaveText* norm_title = new TPaveText(xCenter,yCenter,xCenter,yCenter, "br");
    norm_title->SetBorderSize(0);
    norm_title->SetFillColor(0);
    norm_title->SetFillStyle(0);
    TText* norm_title_text = norm_title->AddText("Normal Ordering");
    norm_title_text->SetTextColor(gNor->GetFillColor());
    norm_title_text->SetTextFont(22);
    norm_title_text->SetTextSize(0.04);

    //norm_title->Draw("same l");
    
    cc->Modified();
    cc->Update();	
  

    ///hFormat->SetTitle(axisTitle.Data());
    
    TGaxis *axis = new TGaxis(cc->GetUxmax(), cc->GetUymax(), cc->GetUxmax()*0.999, cc->GetUymin(), 
                              1000.0*sqrt(A*1.0e24/fYmax), 1000.0*sqrt(A*1.0e24/fYmin), 510, "-LG");
    //axis->SetTitleOffset(1.1);
    axis->SetLabelFont(22);
    axis->SetLabelSize(0.05);
    axis->SetLabelOffset(0.005);
    axis->SetTitleFont(22);
    axis->SetTitleSize(0.06);
    axis->SetTitleOffset(0.85);
    
    axis->CenterTitle();
    axis->SetTitle("m_{#beta#beta} (meV)");
    //axis->SetTitle(axisTitle.Data());

    cc->cd();
    axis->Draw();
  }
  
  cc->SetLogy();
  cc->SetTickx();
  cc->RedrawAxis();
  cc->Modified();
  cc->Update();
  
  return cc;
}

void nEXOHalflifeVsTimePlot::PlotGraph(TGraph& graph, TCanvas& canvas, TLegend& leg, TF1* smooth)
{
  if(smooth)
  {
    smooth->SetRange(0.9*fXmin,1.1*fXmax);
    for(int i = 0; i < 10; i++)
      graph.Fit(smooth,"QR");
    //graph.Fit(smooth,"R");
    smooth->SetLineColor(graph.GetLineColor());
    smooth->SetLineWidth(graph.GetLineWidth());
    smooth->SetLineStyle(graph.GetLineStyle());
    smooth->SetNpx(10000);
    canvas.cd();
    smooth->Draw("same");
    leg.AddEntry(smooth,graph.GetTitle(),"L");
  }
  else
  {
    canvas.cd();
    graph.Draw("same l");
    leg.AddEntry(&graph,graph.GetTitle(),"L");
  }
}


void nEXOHalflifeVsTimePlot::PlotEXO200(TCanvas& canvas, TLegend& leg)
{  
  Double_t hl_EXO200_d50 = 2.6e25;
  Double_t hl_EXO200_sen = 5.7e25;
  
  TLine* lineEXO200_d50 = new TLine(fXmin, hl_EXO200_d50, fXmax, hl_EXO200_d50);//new TLine(0., hl_EXO200_d50, 3.5, hl_EXO200_d50);
  //TArrow* lineEXO200_d50 = new TArrow(0., hl_EXO200_d50, 3., hl_EXO200_d50, 0.01, "<");
  ////TLine* lineEXO200_d90 = new TLine(0., hl_EXO200_d90, 3., hl_EXO200_d90);
  TLine* lineEXO200_sen = new TLine(fXmin, hl_EXO200_sen, fXmax, hl_EXO200_sen);//new TLine(0., hl_EXO200_sen, 3.5, hl_EXO200_sen);
  //TArrow* lineEXO200_sen = new TArrow(0., hl_EXO200_sen, 3., hl_EXO200_sen, 0.01, "<");
  
  lineEXO200_d50->SetLineColor(kRed);
  //lineEXO200_d90->SetLineColor(kMagenta);
  lineEXO200_sen->SetLineColor(kBlue);	
  
  lineEXO200_d50->SetLineWidth(3);
  //lineEXO200_d50->SetLineStyle(10);
  //lineEXO200_d90->SetLineWidth(2);
  lineEXO200_sen->SetLineWidth(3);
  
  lineEXO200_d50->SetLineWidth(3);
  //lineEXO200_d90->SetLineStyle(10);
  lineEXO200_sen->SetLineWidth(3);
  
  lineEXO200_d50->SetLineStyle(7);//kDashed);
  //lineEXO200_d90->SetLineStyle(kDashed);
  lineEXO200_sen->SetLineStyle(7);//kDashed);
  
  TPaveText* exo200Ultimate_title = new TPaveText(0.25, 0.205, 0.25, 0.205, "brNDC");
  exo200Ultimate_title->SetBorderSize(0);
  exo200Ultimate_title->SetFillColor(0);
  exo200Ultimate_title->SetFillStyle(0);
  TText* exo200Ultimate_title_text = exo200Ultimate_title->AddText("EXO200, Ultimate");
  exo200Ultimate_title_text->SetTextColor(kBlack);
  exo200Ultimate_title_text->SetTextFont(22);
  exo200Ultimate_title_text->SetTextSize(0.04);

  canvas.cd();
  lineEXO200_sen->Draw("same l");
  leg.AddEntry(lineEXO200_sen,"EXO-200 Phase-II Sensitivity (90% CL)","L");
  //exo200Ultimate_title->Draw("same l");
  lineEXO200_d50->Draw("same l");
  leg.AddEntry(lineEXO200_d50,"EXO-200 Phase-II Discovery 3#sigma (50% Prob.)","L");
  //lineEXO200_d90->Draw();
}
  
