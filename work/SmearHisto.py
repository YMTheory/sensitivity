
import ROOT
import copy, math
import os,sys
from optparse import OptionParser

verbose = 1

def SmearHisto1D(histoIn, resol, rebin = 1):

    binIn = histoIn.GetXaxis()
    elow = binIn.GetBinLowEdge(1)

    #histoOut = copy.copy(histoIn)
    histoOut = histoIn.Rebin(rebin,histoIn.GetName()+'_GausSmear')
    histoOut.Reset("ICESM")
    binOut = histoOut.GetXaxis()

    for bi in range(1,binIn.GetNbins()+1):
        if verbose >= 2 and bi % 100 == 0:
            print bi, binIn.GetNbins()
        emc = binIn.GetBinCenter(bi)
        res = resol.Eval(emc)
        constres = 0.7071067811865474 / res
        erflow = math.erf(constres * ( elow - emc ) )
        weight = histoIn.GetBinContent(bi)
        for bo in range(1,binOut.GetNbins()+1):
            erfup = math.erf(constres * (binOut.GetBinUpEdge(bo) - emc) )
            histoOut.SetBinContent(bo,histoOut.GetBinContent(bo) + (erfup - erflow) * weight)
            erflow = erfup
    
    histoOut.Scale(0.5)
    
    return histoOut

def SmearHisto2D(histoIn, resol, rebin = 1):
    
    binInY = histoIn.GetYaxis()

    histoOut = histoIn.Rebin2D(rebin,1,histoIn.GetName()+'_Smear')
    histoOut.Reset("ICESM")

    for by in range(1,binInY.GetNbins()+1):
        if verbose >= 1:
            print by, binInY.GetNbins()
        histoIn1D = histoIn.ProjectionX(histoIn.GetName()+'_1d', by, by)
        histoOut1D = SmearHisto1D(histoIn1D, resol, rebin)
        binOutX = histoOut1D.GetXaxis()

        for bx in range(1,binOutX.GetNbins()+1):
            histoOut.SetBinContent(bx,by,histoOut1D.GetBinContent(bx))

    return histoOut  


if __name__ == "__main__":

    libRelPath = '../lib/libnEXOSensitivity.so'
    ROOT.gSystem.Load(libRelPath)

    usage = "usage: python RunSensitivity.py [options]"
    parser = OptionParser(usage)
    
    parser.add_option("-r","--resolution", nargs=1,type=float,default=1.)
    parser.add_option("-i","--infile", nargs=1)
    parser.add_option("-o","--outfile", nargs=1)
    #parser.add_option("-p","--pdf", nargs=1)
    parser.add_option("-b","--rebin", nargs=3,type=int)    

    options,args = parser.parse_args()   
    print 'Using options:', options


    # set resolution function
    p0, p1, p2 = 0., 36.9, 8.8e-3
    q = 2458.

    resol = ROOT.TF1('resol',"[3]*TMath::Sqrt([0]*[0]*x + [1]*[1] + [2]*[2]*x*x)",0,100000)
    resol.SetParameters(p0,p1,p2,1.)
    qResol = resol.Eval(q)/q
    resol.SetParameter(3,options.resolution/qResol)
    
    # get histogram to e smeared
    inFile = ROOT.TFile.Open(options.infile)
    inHistSS = inFile.Get('h_StandoffVsEnergySS')
    inHistMS = inFile.Get('h_StandoffVsEnergyMS')
    
    inHistSS.Rebin2D(options.rebin[0],options.rebin[1])
    inHistMS.Rebin2D(options.rebin[0],options.rebin[1])

    print 'Smearing SS ...'
    outHistSS = SmearHisto2D(inHistSS,resol,options.rebin[2])
    print 'Smearing MS ...'
    outHistMS = SmearHisto2D(inHistMS,resol,options.rebin[2])

    print 'Writing smeared histograms ...'
    outFile = ROOT.TFile.Open(options.outfile,'recreate')
    outHistSS.Write()
    outHistMS.Write()
    outFile.Close()

    print 'Saved into file:', options.outfile

    inFile.Close()

