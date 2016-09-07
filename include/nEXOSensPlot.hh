#ifndef nEXOSensPlot_hh
#define nEXOSensPlot_hh

#include <iostream>

#include "TNamed.h"

class nEXOSensPlot : public TNamed
{
public:
  nEXOSensPlot(const char* name = 0, const char* title = 0);
  virtual ~nEXOSensPlot();

  virtual TObject* GetPlot() = 0;
  
protected:
  
  ClassDef(nEXOSensPlot,1) // Base class for nEXO sensitivity plots
};

#endif

