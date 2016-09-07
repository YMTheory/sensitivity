#ifndef nEXONuclearMatrixElement_hh
#define nEXONuclearMatrixElement_hh

#include <iostream>
#include <vector>

#include "TString.h"
#include "TNamed.h"

class nEXONuclearMatrixElement : public TNamed
{
public:
  enum NME_t {kGCM=0, kQRPA2, kRQRPA, kNSM, kIBM2};
  static const size_t nNMEs = 5;
  
  nEXONuclearMatrixElement(NME_t nme);
  virtual ~nEXONuclearMatrixElement();

  Double_t GetAxisScale(){return fAxisScale;}
  const char* GetAxisTitle(){return fAxisTitle.Data();}

private:

  const NME_t fNME;

protected:

  Double_t fAxisScale;
  TString fAxisTitle;

  ClassDef(nEXONuclearMatrixElement,1)	// Class for NME
};

#endif
