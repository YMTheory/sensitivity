

import ROOT
import sys

libRelPath = '../lib/libnEXOSensitivity.so'
ROOT.gSystem.Load(libRelPath)

print ROOT.nEXOSensitivity.kRdmCV,  int(ROOT.nEXOSensitivity.kUL) #kPosUL #kRdmCV
