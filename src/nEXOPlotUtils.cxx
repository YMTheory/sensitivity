//////////////////////////////////////////////////////////////////////////////////
//
//
// nEXOPlotUtils is a namespace with a collection of functions useful for plotting
//
//
//////////////////////////////////////////////////////////////////////////////////

#include "nEXOUtils.hh"

/*
TCanvas* nEXOPlotUtils::CreateHalflifeVsTimeCanvas(nEXONuclearMatrixElement::NME_t matElem, const char* name, Double_t xMin, Double_t xMax, Double_t yMin, Double_t yMax)
{
  nEXONuclearMatrixElement nme(matElem);
  
  TString axisTitle = nme.GetTitle();
  Double_t A = nme.GetAxisScale();

  TString canvasName(name);
  canvasName.ReplaceAll("[NME]",nme.GetName());

  TString hFormatName = Form("hForm_%s", nme.GetName());
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
       
  return 0;
}
*/
