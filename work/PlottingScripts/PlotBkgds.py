
import ROOT

ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')

filename = '../results/done/fits_db_v73_2016-09-26_bb2n_0nu_rdm_10.0_years_0.0_counts_*.root' # '../results/done/fits_hamamatsu_v68_2016-06-21_0nu_red1x_fine_rdm_5.0_years_0.0_counts_*.root'

def GetBkgdLimit(chain,branch,prob,max = 0.10):

    chain.SetEstimate(chain.GetEntries()+1)
    branch = branch + '>>hist(%d,0,%d)'%(max*1000,max)
    chain.Draw(branch,'','goff')
    #print 'using', chain.GetSelectedRows()
    hist = ROOT.gDirectory.Get('hist')
    hist.Scale(1./hist.Integral())
    probs = ROOT.TGraph()
    for b in range(1,hist.GetNbinsX()+1):
        probs.SetPoint(probs.GetN(),hist.Integral(1,b),hist.GetBinCenter(b))
    
    return probs.Eval(prob)


chain = ROOT.TChain('tree')
chain.Add(filename)

canvas = ROOT.TCanvas()
canvas.Divide(3,1)

canvas.cd(1)
chain.Draw('bkg_fwhm_1t>>h1(100,0,0.25)','','goff')
h1 = ROOT.gDirectory.Get('h1')
h1.Scale(1./h1.Integral())
h1.GetXaxis().SetTitle('Background in FWHM-1t (counts/year)')
h1.GetYaxis().SetTitle('Toy Fits (normalized)')
h1.Draw()
canvas.cd(2)
chain.Draw('bkg_fwhm_3t>>h3(100,0,5)','','goff')
h3 = ROOT.gDirectory.Get('h3')
h3.Scale(1./h3.Integral())
h3.GetXaxis().SetTitle('Background in FWHM-3t (counts/year)')
h3.GetYaxis().SetTitle('Toy Fits (normalized)')
h3.Draw()
canvas.cd(3)
chain.Draw('bkg_fwhm_fv>>h(100,0,10)','','goff')
h = ROOT.gDirectory.Get('h')
h.Scale(1./h.Integral())
h.GetXaxis().SetTitle('Background in FWHM-FV (counts/year)')
h.GetYaxis().SetTitle('Toy Fits (normalized)')
h.Draw()

print GetBkgdLimit(chain,'bkg_fwhm_1t',0.95), GetBkgdLimit(chain,'bkg_fwhm_3t',0.95), GetBkgdLimit(chain,'bkg_fwhm_fv',0.95)
print GetBkgdLimit(chain,'bkg_fwhm_1t',0.5), GetBkgdLimit(chain,'bkg_fwhm_3t',0.5), GetBkgdLimit(chain,'bkg_fwhm_fv',0.5)
print GetBkgdLimit(chain,'bkg_fwhm_1t',0.9), GetBkgdLimit(chain,'bkg_fwhm_3t',0.9), GetBkgdLimit(chain,'bkg_fwhm_fv',0.9)

#chi3 = ROOT.TF1('chi3','[0]*ROOT::Math::chisquared_pdf(x,[1])',0,10)
#chi3.SetParameters(10,10)
#h3.Fit(chi3)

raw_input('')
