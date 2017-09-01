#include <iostream>

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH2.h"
#include "TMatrixDSym.h"

using namespace RooFit;



double fitUL = 23.92;

double getMedian(TH1D* h1, double& lo, double& hi) {
    double xq[]={.5,0.025,0.975};
    double yq[]={.5,0.025,0.975};
    h1->GetQuantiles(3,xq,yq);
    double median = xq[0];
    lo = xq[1]; hi = xq[2];
    return median;
}


void PlotBb0nUL_ts(Int_t withBb0n=false) {
  TString inFileName = "/p/lscratche/nexouser/sens_recalc/results/done/*.root"; //path to the file containing the fit results
    TChain ch("tree");
    ch.Add(inFileName.Data());
   
    TH1D* h_ul  = new TH1D("h_up","90% UL on #beta#beta0#nu Counts", 1000, 0., 50.);        

    h_ul->GetXaxis()->SetTitle("Upper Limit on #beta#beta0#nu Counts (90% CL)");

    //set the addresses of the wanted variables
    ch.Draw("(num_signal + num_signal_eLo)>>h_up","nll_ratio>=-.01 && stat_sig==0 && stat_bkg==0 && covQual_sig==3 && covQual_bkg==3","goff");
    // Double_t num_bb0n_eHi = 0.;
    // ch.SetBranchStatus("*", 0);
    // ch.SetBranchStatus("num_signal_eHi", 1);
    
    // ch.SetBranchAddress("num_signal_eHi", &num_bb0n_eHi);
    // //ch.SetBranchAddress("num_bb0n_ul", &num_bb0n_ul);

    // //ch.SetBranchAddress("fitres_sig",&fitres);

    // //Fill the histograms
    // for (int i=0; i<ch.GetEntries(); i++) {
    //       //cout << i << " " << ch.GetEntries() << endl;
    //     ch.GetEntry(i);

    //     h_ul->Fill(num_bb0n_eHi);  
        
    // }

    //construct 2D color plot of correlation matrix
    //gStyle->SetOptStat(0);
    //gSTyle->SetPalette(1);
    //TH2* hcorr = ch.CorrelationHist();
    
    // Extract correlation matrix as TMatrixDSym
    //const TMatrixDSym& cor = ch.correlationMatrix();
    //Print correlation matrix
    //cout<<"correlation matrix" << endl;
    //cor.Print();

    //calculate the mean, median, and 68% interval containing the median
    Double_t loE = 0., hiE = 0.;
    Double_t mean = h_ul->GetMean();
    Double_t median = getMedian(h_ul, loE, hiE);
    cout << "Mean counts = " << mean << ", median 90% UL = " << median << endl;
    cout << loE << " , " << hiE << endl;
    

    //Find the bins for the 68% interval
    h_ul->Rebin(5);
    Int_t loBin = h_ul->GetXaxis()->FindBin(loE);
    Int_t hiBin = h_ul->GetXaxis()->FindBin(hiE);

    //Histo to display the 68% interval around the median
    TH1D* h_interval = new TH1D("h_interval", ";Upper Limit on #beta#beta0#nu Counts (90% CL);Number of Toy Experiments", 200, 0., 50.);
    for (int i=loBin; i<=hiBin; i++) {
        h_interval->SetBinContent(i, h_ul->GetBinContent(i));
    }
    h_interval->SetFillColor(kCyan);
    h_interval->SetFillStyle(3001);
    h_interval->SetLineColor(0);
    h_interval->GetYaxis()->SetRangeUser(0,2500);
    h_interval->GetYaxis()->SetTitleOffset(1.25);
    //Find the bins corresponding to the median upper limit, and the upper limit from real data
    Double_t medBinCon = h_ul->GetBinContent(h_ul->GetXaxis()->FindBin(median));
    Double_t fitBinCon = h_ul->GetBinContent(h_ul->GetXaxis()->FindBin(fitUL)); 

    //Draw the results
    TCanvas* cc_num = new TCanvas("cc","",1200,900); 
    h_interval->DrawCopy("")->SetStats(0);
    h_ul->DrawCopy("same");
    TLine* medianLine = new TLine(median, 0., median, medBinCon);
    medianLine->SetLineColor(kRed);
    medianLine->SetLineWidth(2);
    //TLine* fitMedianLine = new TLine(fitUL, 0., fitUL, fitBinCon);
    //fitMedianLine->SetLineColor(kBlue);
    //fitMedianLine->SetLineWidth(2);

    TLegend* leg = new TLegend(0.475, 0.45, 0.9, 0.7);
    //leg->AddEntry(fitMedianLine, "Upper limit from data, 23.92 counts", "l");
    TString entryTitle = "Median upper limit"; entryTitle += Form(", %0.2f counts", median);
    leg->AddEntry(medianLine, entryTitle, "l");
    leg->AddEntry(h_interval, "2#sigma Interval around Median", "f");

    medianLine->Draw();
    //fitMedianLine->Draw();
    leg->Draw();
    gPad->RedrawAxis();
    cc_num->SaveAs("/p/lscratche/nexouser/sens_recalc/work/countlimits.pdf");
    cc_num->SaveAs("/p/lscratche/nexouser/sens_recalc/work/countlimits.root");
}
