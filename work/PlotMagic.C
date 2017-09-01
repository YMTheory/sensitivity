Double_t spline_func(Double_t *x, Double_t *par)
{
    /*Fit parameters:
     par[0-3]=X of nodes (to be fixed in the fit!)
     par[4-7]=Y of nodes
     par[8-9]=first derivative at begin and end (to be fixed in the fit!)
     */
    Double_t xx = x[0];
    const Int_t nknots = 9;
    
//    Double_t xn[5] = { par[0], par[1], par[2], par[3], par[4] };
//    Double_t yn[5] = { par[5], par[6], par[7], par[8], par[9] };
    Double_t xn[nknots] = { 0.1, 1., 5., 7., 10., 15, 20, 30, 50. };
    Double_t yn[nknots] = { par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8] };
    
//    Double_t b1 = par[8];
//    Double_t e1 = par[9];
    
//    TSpline3 sp3("sp3", xn, yn, 4, "b1e1", b1, e1);
    TSpline3 sp3("sp3", xn, yn, nknots, "b1e1", 0, 0);
    
    return sp3.Eval(xx);
}


void PlotMagic()
{
    
    string basename = "../results/all_bkgs-1DoF/allbkgs_10.0_years_";
    //std::vector<Double_t> numcounts_arr = {0.1,2,5,7,9,9.1,9.2,9.3,9.4,9.5,9.6,9.7,9.8,9.9,10,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9,11,11.1,11.2,11.3,11.4,11.5,11.6,11.7,11.8,11.9,12,12.1,12.2,12.3,12.4,12.5,12.6,12.7,12.8,12.9,13,13.1,13.2,13.3,13.4,13.5,13.6,13.7,13.8,13.9,14,14.1,14.2,14.3,14.4,14.5,14.6,14.7,14.8,14.9,15,20,25,30};
    std::vector<Double_t> numcounts_arr = {0.1, 1, 2};
    std::vector<Double_t>  magic_numbers;
    
    for (int i=5; i <= 30; i++) numcounts_arr.push_back(i);
    numcounts_arr.push_back(50);
    
    for (int i = 0; i < numcounts_arr.size(); ++i) {
        
        Double_t numcounts = numcounts_arr[i];
        
        TChain *chain = new TChain("tree");
        chain->Add(Form("%s%.2f_*", basename.c_str(), numcounts));
        
        TH1F ratiohist("ratiohist",Form("ratiohist_%.2f",numcounts), 800, 0, 20);
        chain->Draw("nll_ratio>>ratiohist",
			"nll_ratio>=0 && stat_sig==0 && stat_bkg==0 && covQual_sig==3 && covQual_bkg==3","goff");
        
        Double_t xq[] = {0.9};  // position where to compute the quantiles in [0,1]
        Double_t yq[] = {0};  // array to contain the quantiles
        ratiohist.GetQuantiles(1,yq,xq);
        
        magic_numbers.push_back(yq[0]/2.);
        cout << numcounts << " " << magic_numbers.back() << " " << ratiohist.GetEntries() << endl;
            
    }
    
    
    auto xmin = 0.;
    auto xmax = 50.;
    Int_t npars = 9;
    
    TF1 *f_spline4 = new TF1("f_spline", spline_func, xmin, xmax, npars); // npars = 2*nodes
    
//    f_spline4->SetParameters(5., 10., 15., 20., 1., 1., 1., 1.);
//    f_spline4->FixParameter(0,numcounts_arr.front());
//    f_spline4->SetParLimits(1,xmin,xmax);
//    f_spline4->SetParLimits(2,xmin,xmax);
//    f_spline4->SetParLimits(3,xmin,xmax);
//    f_spline4->FixParameter(4,numcounts_arr.back());
//    f_spline4->FixParameter(5,magic_numbers.front());
//    f_spline4->SetParLimits(6,1,3);
//    f_spline4->SetParLimits(7,1,3);
//    f_spline4->SetParLimits(8,1,3);
//    f_spline4->FixParameter(9,magic_numbers.back());

    f_spline4->SetParameters(2,2,2,2,2,2,2,2,2);
//    f_spline4->FixParameter(0,magic_numbers.front());
//    f_spline4->FixParameter(8,magic_numbers.back());
    
    
    TGraph *gr = new TGraph(magic_numbers.size(), numcounts_arr.data(), magic_numbers.data());
    
    gr->Draw("APL");
    
    gr->Fit(f_spline4, "R");
    
    cout << f_spline4->Eval(0) << endl;

}
