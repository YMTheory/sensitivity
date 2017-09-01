
import ROOT
import copy, math
import SmearHisto as sh

# verbose
sh.verbose = 0

# input values
q_value = 2458
e_thres = 1000
bb2n_decays_3t_ss = 2945466.3 *3/3.74 # bb2n decays in inner 3 tonne = total * 3/3.74
bb2n_decays_3t_ms = 42845.0*3/3.74
bb2n_tot_yr = bb2n_decays_3t_ss #2406782.654
bb2n_roi_eff = 1 #0.813
sd_cut = 9 # cm
vol = {9: 3}  # in tonne
resqs = [0.01041] #0.005,0.01,0.015,0.020,0.025,0.029,0.047]

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
# hin.Add(fin.Get('h_StandoffVsEnergyMS'))
# bb2n_tot_yr += bb2n_decays_3t_ms
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

# results for SS-only without efficiency correction
# nexo {0.005: {'1-sigma': 9.05239438466765e-05, '3-sigma': 0.0025381586041389598, 'fwhm': 0.00012607273026077337, 'klz': 4.749505188418214, '2-sigma': 0.000561590648584076}, 0.025: {'1-sigma': 1.042865445080679, '3-sigma': 30.56405429411815, 'fwhm': 1.5021475556859514, 'klz': 16.31195528459692, '2-sigma': 6.631785862716242}, 0.015: {'1-sigma': 0.05342835932060552, '3-sigma': 1.6396064587808343, 'fwhm': 0.07885523332515731, 'klz': 7.99912023002309, '2-sigma': 0.3518082155999309}, 0.01: {'1-sigma': 0.005075871587905567, '3-sigma': 0.15323328465320318, 'fwhm': 0.007090731185826371, 'klz': 5.895668336819473, '2-sigma': 0.03417090191979355}, 0.047: {'1-sigma': 35.58736350480467, '3-sigma': 978.8890842966052, 'fwhm': 51.33279811707325, 'klz': 72.91405688742816, '2-sigma': 224.65753473762015}, 0.029: {'1-sigma': 2.4363004183396697, '3-sigma': 69.47918140400645, 'fwhm': 3.4317764120351057, 'klz': 21.801958044036994, '2-sigma': 15.473860577332744}, 0.02: {'1-sigma': 0.2933612628476112, '3-sigma': 8.609958657755293, 'fwhm': 0.4091890067429631, 'klz': 11.34022394771196, '2-sigma': 1.8867066754005464}}
# klzu {0.005: {'1-sigma': 8.916245175782933e-05, '3-sigma': 0.0025118795949373163, 'fwhm': 0.00012423609560840987, 'klz': 4.734572211058359, '2-sigma': 0.0005546056412137901}, 0.025: {'1-sigma': 0.9704655996174552, '3-sigma': 29.157820448232684, 'fwhm': 1.4015645420877263, 'klz': 15.492513683062393, '2-sigma': 6.257088229439432}, 0.015: {'1-sigma': 0.05109944418290979, '3-sigma': 1.5916061510875625, 'fwhm': 0.07554383644219342, 'klz': 7.813562390284531, '2-sigma': 0.3392699151548868}, 0.01: {'1-sigma': 0.004925656455498029, '3-sigma': 0.15015938695428144, 'fwhm': 0.006887513413175839, 'klz': 5.82751864584081, '2-sigma': 0.03334517130667791}, 0.047: {'1-sigma': 31.331163711845875, '3-sigma': 903.2896123427306, 'fwhm': 45.41658625146374, 'klz': 64.8201665654633, '2-sigma': 202.96758445951855}, 0.029: {'1-sigma': 2.2437309461529367, '3-sigma': 65.83123434646916, 'fwhm': 3.1696390034630895, 'klz': 20.46099674711195, '2-sigma': 14.479705705525703}, 0.02: {'1-sigma': 0.27672085314407013, '3-sigma': 8.284039526274668, 'fwhm': 0.38671753671587794, 'klz': 10.927773644192776, '2-sigma': 1.7993643707266074}}

