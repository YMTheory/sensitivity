#include "nEXONuclearMatrixElement.hh"

ClassImp(nEXONuclearMatrixElement);

nEXONuclearMatrixElement::nEXONuclearMatrixElement(NME_t nme) : TNamed(), fNME(nme)
{

  // set name, title and scale
  switch(fNME)
  {
      case kGCM:
        SetName("gcm");
        SetTitle("GCM");
        fAxisScale = 0.31;
        break;
      case kQRPA2:
        SetName("qrpa2");
        SetTitle("QRPA-2");
        fAxisScale = 2.21;
        break;
      case kRQRPA:
        SetName("rqrpa");
        SetTitle("R-QRPA");
        fAxisScale = 0.83;
        break;
      case kNSM:
        SetName("nsm");
        SetTitle("NSM");
        fAxisScale = 1.78;
        break;
      case kIBM2:
        SetName("ibm2");
        SetTitle("IBM-2");
        fAxisScale = 0.49;
        break;
  }

  fAxisTitle = Form("m_{#beta#beta} (meV) %s",GetName());
}

nEXONuclearMatrixElement::~nEXONuclearMatrixElement()
{
}
