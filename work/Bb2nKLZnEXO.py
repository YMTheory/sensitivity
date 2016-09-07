
import ROOT
import copy, math
import SmearHisto as sh

# verbose
sh.verbose = 0

# input values
q_value = 2458
e_thres = 1000
bb2n_tot_yr = 2406782.654
bb2n_roi_eff = 0.813
sd_cut = 9 # cm
vol = {9: 3}  # in tonne
resqs = [0.005,0.01,0.015,0.020,0.025,0.029,0.047]

# define roi and limits
rois = {'1-sigma':1,'2-sigma':2,'3-sigma':3,'fwhm':(2*math.log(2))**0.5,'klz':0}
roilos = {}
roihis = {}
roilos['klz'] = 2300
roihis['klz'] = 2700

# resolution functions
nexo_resol = ROOT.TF1('nexo_resol',"[3]*TMath::Sqrt([0]*[0]*x + [1]*[1] + [2]*[2]*x*x)",0,100000)
nexo_resol.SetParameters(0.,36.9,8.8e-3,1./1.74)
klz_resol = ROOT.TF1('klz_resol','[0]/TMath::Sqrt(x/1000.)*x',0,10000)
klz_resol.SetParameter(0,0.045)
print 'resol nexo & klz:', nexo_resol.Eval(2458), klz_resol.Eval(2458)

# open file and get unsmeared histo
fin = ROOT.TFile.Open('../histos/bb2n/nEXO_Histos_bb2n.root')
hin = copy.copy(fin.Get('h_StandoffVsEnergySS'))
hin.Add(fin.Get('h_StandoffVsEnergyMS'))
hene = hin.ProjectionX('%s_ene'%(hin.GetName()),sd_cut,-1)

# evaluate threshold normalization based on unsmeared histogram
thres_norm = hene.Integral(hene.FindBin(e_thres),-1)/hene.Integral()
print 'threshold normalization', thres_norm

# auxiliary function to change energy threshold to get rid of edge effect when smearing
def ChangeThreshold(histo,thres):

    bt = histo.FindBin(thres)
    nb = histo.GetNbinsX() - bt + 1
    nh = ROOT.TH1F(histo.GetName(),histo.GetTitle(),nb,thres,histo.GetXaxis().GetXmax())

    for b in range(bt,histo.GetNbinsX()+1):
        nh.SetBinContent(b-bt+1,histo.GetBinContent(b))

    return nh

# smear histograms for all resolutions
smeared_file = ROOT.TFile.Open('klzu_vs_nexo_bb2n.root','recreate')
nexo_results = {}
klzu_results = {}
for resq in resqs:
    print 'working on resolution', resq
    nexo_results[resq] = {}
    klzu_results[resq] = {}
    
    # adjust resolution functions
    nexo_par = resq/0.0174
    klz_par = resq * (q_value/1000.)**0.5
    nexo_resol.SetParameter(3,nexo_par)
    klz_resol.SetParameter(0,klz_par)

    # smear histograms
    hnexo = sh.SmearHisto1D(hene, nexo_resol, 1)
    hklzu = sh.SmearHisto1D(hene, klz_resol, 1)

    # increase thresholds
    hnexo = ChangeThreshold(hnexo,e_thres)
    hklzu = ChangeThreshold(hklzu,e_thres)

    # normalize to expected bb2n counts per tonne per year
    hnexo.Scale(1./hnexo.Integral() * thres_norm * bb2n_tot_yr / bb2n_roi_eff / vol[sd_cut])
    hklzu.Scale(1./hklzu.Integral() * thres_norm * bb2n_tot_yr / bb2n_roi_eff / vol[sd_cut])

    smeared_file.cd()
    hnexo.SetName('%s_nexo_res%.3f'%(hnexo.GetName(),resq))
    hnexo.Write()
    hklzu.SetName('%s_klzu_res%.3f'%(hklzu.GetName(),resq))
    hklzu.Write()

    # evaluate bb2n for each roi
    sigma = resq*q_value
    for roi in rois:
        if not roi in ['klz']:
            roilos[roi] = q_value - sigma*rois[roi]
            roihis[roi] = q_value + sigma*rois[roi]  

    for roi in rois:
        print 'evaluating roi', roi, roilos[roi], roihis[roi]
        nexo_res = hnexo.Integral(hnexo.FindBin(roilos[roi]),hnexo.FindBin(roihis[roi]))
        klzu_res = hklzu.Integral(hklzu.FindBin(roilos[roi]),hklzu.FindBin(roihis[roi]))
        print 'nexo', nexo_res, 'klz', klzu_res
        nexo_results[resq][roi] = nexo_res
        klzu_results[resq][roi] = klzu_res

print 'nexo', nexo_results
print 'klzu', klzu_results    

fin.Close()

# results
# nexo {0.005: {'1-sigma': 0.00012156702423737897, '3-sigma': 0.0033944465261427914, 'fwhm': 0.0001692518237632612, 'klz': 6.297484902155599, '2-sigma': 0.0007526356585967342}, 0.025: {'1-sigma': 1.3823112294776365, '3-sigma': 40.42842561928248, 'fwhm': 1.990805720575736, 'klz': 21.588934716815682, '2-sigma': 8.782721038129239}, 0.015: {'1-sigma': 0.07094488432994694, '3-sigma': 2.174251918499273, 'fwhm': 0.10469339584233239, 'klz': 10.598819211820723, '2-sigma': 0.4668243877624718}, 0.01: {'1-sigma': 0.0067560696852524416, '3-sigma': 0.20346954650063775, 'fwhm': 0.009434754255835287, 'klz': 7.8150010856290635, '2-sigma': 0.04540726057558686}, 0.047: {'1-sigma': 46.97428731061518, '3-sigma': 1285.7508966185196, 'fwhm': 67.73201160365716, 'klz': 96.1696785929671, '2-sigma': 295.85924841600354}, 0.029: {'1-sigma': 3.227244202804286, '3-sigma': 91.78470542957535, 'fwhm': 4.545140399015509, 'klz': 28.837731118066934, '2-sigma': 20.473612865109317}, 0.02: {'1-sigma': 0.3891544728685403, '3-sigma': 11.404611063743879, 'fwhm': 0.5427439321938436, 'klz': 15.018502268045767, '2-sigma': 2.5010184130460402}}
# klzu {0.005: {'1-sigma': 0.00011974399735947827, '3-sigma': 0.0033594812077865424, 'fwhm': 0.00016679324713209098, 'klz': 6.277691980182578, '2-sigma': 0.0007433006020418631}, 0.025: {'1-sigma': 1.2863696686399635, '3-sigma': 38.57113451546681, 'fwhm': 1.8575403469003504, 'klz': 20.505626336916766, '2-sigma': 8.286820582914515}, 0.015: {'1-sigma': 0.06785368563942029, '3-sigma': 2.1106198985184577, 'fwhm': 0.10029860512258892, 'klz': 10.35305068180718, '2-sigma': 0.4501921377859617}, 0.01: {'1-sigma': 0.006556537740834756, '3-sigma': 0.19938936392624984, 'fwhm': 0.009164839125332946, 'klz': 7.724708317088643, '2-sigma': 0.04431071400827591}, 0.047: {'1-sigma': 41.364399360027164, '3-sigma': 1186.7464309277852, 'fwhm': 59.93791936291382, 'klz': 85.51202678386471, '2-sigma': 267.35458876205666}, 0.029: {'1-sigma': 2.972276038664859, '3-sigma': 86.97410939495074, 'fwhm': 4.198150982731022, 'klz': 27.066324976875123, '2-sigma': 19.1597393540942}, 0.02: {'1-sigma': 0.36708573205396533, '3-sigma': 10.973135812079157, 'fwhm': 0.5129452611581655, 'klz': 14.47264017687168, '2-sigma': 2.3852696189496783}}
