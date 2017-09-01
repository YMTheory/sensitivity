
import ROOT
import sys, array

#nmes = [ROOT.nEXONuclearMatrixElement.kGCM, ROOT.nEXONuclearMatrixElement.kQRPA2, ROOT.nEXONuclearMatrixElement.kRQRPA, ROOT.nEXONuclearMatrixElement.kNSM, ROOT.nEXONuclearMatrixElement.kIBM2] 

def PlotMbb():

    rout = ROOT.TFile('mbb_vs_mnu.root','recreate')

    for nmo in range(ROOT.nEXOMbbVsMnuPlot.nNuMassOptions):
        plot = ROOT.nEXOMbbVsMnuPlot('cc')
        plot.fNuMassOption = nmo
        plot.SetBandFromFile("nexo_sens_gamma","nEXO, No Ba-tagging",'../results/full/allfits_hamamatsu_0nu_eff_tpc_rdm_5.0_years_0.0_counts.root',5,3740)
        #plot.SetBand("nexo_sens_gamma","nEXO, No Ba-tagging",0.0072, 0.0192)
        #plot.SetBand("nexo_sens_batag","nEXO, With Ba-tagging",0.0031, 0.0083)
        ROOT.nEXONuOscPars.GetInstance().SetFitSource(ROOT.nEXONuOscPars.kprd90_093006_2014)
        plot.fLoXRangeCoeff = 5
        plot.fLoXRangeExp = 2
        plot.fLoYRangeExp = 3
        canvas = plot.GetPlot()
        
        rout.cd()
        canvas.Write()

    rout.Close()

def PlotHL():

    rout = ROOT.TFile('hlv7.root','recreate')

    yrs = array.array('d',[0.5, 1., 2.5, 5., 10.])#,12.5,15.0,20.0])
    cts = array.array('d',[0,2.5,5,10,20,40,80]) #range(11))
    #count_sens = ROOT.nEXOUtils.EvalCountSensFromFiles(len(yrs),yrs,'../../v5/results/done/disc_hamamatsu_v62_0nu_eff_tpc_rdm_%0.1f_years_0.0_counts_*.root','bkg_fwhm_3t',3000,0.81) # 80% eff to find 0nu in FV > 700 keV, 70% eff to be in FWHM
    #count_disc = ROOT.nEXOUtils.EvalCountDiscFromFiles(len(yrs),yrs,'../../v5/results/done/disc_hamamatsu_v62_0nu_eff_tpc_rdm_%0.1f_years_0.0_counts_*.root','bkg_fwhm_3t',3000,0.81)
    disc = ROOT.nEXOUtils.ReadDiscFromFiles(len(yrs),yrs,len(cts),cts,'../../v5/results/done/disc_hamamatsu_v62_0nu_eff_tpc_rdm_%0.1f_years_%0.1f_counts_*.root',3740,0.5,3)
    sens = ROOT.nEXOUtils.ReadSensFromFiles(len(yrs),yrs,'../../v5/results/done/disc_hamamatsu_v62_0nu_eff_tpc_rdm_%0.1f_years_0.0_counts_*.root',3740)
    #count_sens = ROOT.nEXOUtils.EvalCountSensFromFiles(len(yrs),yrs,'../../v3/results/done/disc_hamamatsu_0nu_eff_tpc_rdm_%0.1f_years_0.0_counts_*.root','bkg_fwhm_3t',3000,0.81) # 80% eff to find 0nu in FV > 700 keV, 70% eff to be in FWHM
    #count_disc = ROOT.nEXOUtils.EvalCountDiscFromFiles(len(yrs),yrs,'../../v3/results/done/disc_hamamatsu_0nu_eff_tpc_rdm_%0.1f_years_0.0_counts_*.root','bkg_fwhm_3t',3000,0.81)
    #disc = ROOT.nEXOUtils.ReadDiscFromFiles(len(yrs),yrs,len(cts),cts,'../../v3/results/done/disc_hamamatsu_0nu_eff_tpc_rdm_%0.1f_years_%0.1f_counts_*.root',3740,0.5,3)
    #sens = ROOT.nEXOUtils.ReadSensFromFiles(len(yrs),yrs,'../../v3/results/done/disc_hamamatsu_0nu_eff_tpc_rdm_%0.1f_years_0.0_counts_*.root',3740)

    for nme in range(ROOT.nEXONuclearMatrixElement.nNMEs):

        print 'Working on NME', nme

        plot = ROOT.nEXOHalflifeVsTimePlot('cc')
        plot.SetNME(nme)
        plot.SetAxisLimits(0,10)

        #yrs = array.array('d',[0.5, 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.])
        #hls = array.array('d',[23.05, 36.04, 56.20, 72.81, 87.44, 100.77, 113.15, 124.77, 135.80, 146.32, 156.41])
        #d50s = array.array('d',[12.55, 19.10, 29.01, 37.02, 44.00, 50.29, 56.10, 61.53, 66.64, 71.51, 76.16])
        #plot.SetGraphPoints('sens_90','Sensitivity (90% CL), 1/10 Cu Bkg',len(yrs),yrs,hls)
        #plot.SetGraphPoints('disc3_50','3#sigma Discovery, 50% Probability, 1/10 Cu Bkg',len(yrs),yrs,d50s)
        #plot.ReadGraphFromFiles('sens_90','Sensitivity (90% CL)',len(yrs),yrs,'../results/done/disc_hamamatsu_0nu_eff_tpc_rdm_%0.1f_years_0.0_counts_*.root',3740)

        plot.SetGraphPoints('sens_90','Sensitivity (90% CL)',len(yrs),yrs,sens,1)
        #plot.SetGraphPoints('count_sens_90','Counting Sensitivity (90% CL), FWHM-3t',len(yrs),yrs,count_sens,1)
        plot.SetGraphPoints('disc3_50','3#sigma Discovery (50% Prob.)',len(yrs),yrs,disc,1)
        #plot.SetGraphPoints('count_disc3_50','Counting 3#sigma Discovery (50% Prob.), FWHM-3t',len(yrs),yrs,count_disc,1)
        

        canvas = plot.GetPlot()

        rout.cd()
        canvas.Write()

    rout.Close()

if __name__ == "__main__":

    ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')

    #PlotMbb()
    PlotHL()
