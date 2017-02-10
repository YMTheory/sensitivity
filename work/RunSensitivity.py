#sample usage: python RunSensitivity.py -n 1000 -r 1 -y 10.0 -c 0.01 -s 1 -b -d ../results -o allbkgs -t ../tables/Summary_v73_2016-09-09_0nu_allbkgs_llnl.root -m 3 --ssfrac-improvement 1.0 --rn222-rate-correction 1.00
#or, for a fine spaced scan: 
#for count in `seq 9 .1 11.9`; do python RunSensitivity.py -n 1000 -r 1 -y 10.0 -c ${count} -s 1 -b -d ../results -o allbkgs -t ../tables/Summary_v73_2016-09-09_0nu_allbkgs_llnl.root -m 3 --ssfrac-improvement 1.0 --rn222-rate-correction 1.00; done
import ROOT
import os,sys
from optparse import OptionParser
import array
if __name__ == "__main__":

    libRelPath = '../lib/libnEXOSensitivity.so'
    ROOT.gSystem.Load(libRelPath)
    # ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

    usage = "usage: python RunSensitivity.py [options]"
    parser = OptionParser(usage)
    
    parser.add_option("--binning", nargs=6,type=float) # choose energy and SD binning, format: n-energy-bins, energy-min, energy-max, n-sd-bins, sd-min, sd-max

    parser.add_option("--ba-tag", action="store_true",default=False) # run in Ba-tag scenario (only bb2n background)
    parser.add_option("--scale-bkgds", nargs=2,type=float) # scale all background counts, 2nd argument whether scale bb2n or not

    parser.add_option("-n","--number-runs", nargs=1,type=int,default=1) # number of toy fits to be produced
    parser.add_option("-s","--random-seed", nargs=1,type=int,default=1) # seed for randomizations
    parser.add_option("-y","--years", nargs=1,type=float,default=5.) # detector livetime in years
    parser.add_option("-t","--tree", nargs=1) # ROOT file with tree for input counts expectation (produced from summary table)
    parser.add_option("-d","--output-dir", nargs=1) # name of the output directory
    parser.add_option("-o","--output-name", nargs=1) # name of the output file

    parser.add_option("-c","--signal-counts", nargs=1,type=float,default=0.0) # included signal counts (for discovery potential)
    parser.add_option("-b","--run-bkgd-only-fit", action="store_true",default=False) # run background only fit (for observation significance)

    parser.add_option("-r","--random-rate", nargs=1,type=int,default=0) # rate (in evt^-1) at which the activity is randomized
    parser.add_option("-m","--method-bkgd", nargs=1,type=int, default=ROOT.nEXOSensitivity.kRdmCV) # method to choose activities : ROOT.nEXOSensitivity.kRdmCV, kUL, kPosUL

    parser.add_option("-g","--group", nargs=1) # modify counts in group given by name here
    parser.add_option("-v","--bkgd-counts", nargs=1,type=float) # modify to counts given here for group given by name above 
    parser.add_option("--ssfrac-improvement", nargs=1,type=float,default=1.) # change of SS/MS fraction for gamma-like PDFs
    parser.add_option("--rn222-rate-correction", nargs=1,type=float,default=1.) # change of Rn222 rate 

    parser.add_option("--turn-off-groups", nargs=1,type=str,default='') # turn these groups off the fit (set zero expectation counts and do not include in fit)
    
    options,args = parser.parse_args()   
    print 'Using options:', options

    realPath = os.path.realpath(options.output_dir)
    realPathWork = realPath + '/working'
    realPathDone = realPath + '/done'
    os.system('mkdir -p %s' % (realPathWork))
    os.system('mkdir -p %s' % (realPathDone))
    outFileName = '%s/%s_%0.1f_years_%0.2f_counts_%04i.root' % (realPathWork,options.output_name,options.years,options.signal_counts,options.random_seed)

    sens = ROOT.nEXOSensitivity(options.random_seed,options.tree)
    sens.fResultFileName = outFileName

    sens.fRunTruthValFit = options.run_bkgd_only_fit
    sens.fExpectCountMethod = options.method_bkgd #ROOT.nEXOSensitivity.kRdmCV #kUL #kPosUL #kRdmCV

    if options.group:
        sens.AddUserMeanCounts(options.group,options.bkgd_counts)

    sens.fVerboseLevel = 0
    sens.fSSFracImprovement = options.ssfrac_improvement
    sens.fRn222RateCorrection = options.rn222_rate_correction

    if options.binning:
        sens.SetBinning(int(options.binning[0]),float(options.binning[1]),float(options.binning[2]),int(options.binning[3]),float(options.binning[4]),float(options.binning[5]))

    sens.SetBaTag(options.ba_tag)
    if options.scale_bkgds:
        sens.ScaleMeanBackgrounds(options.scale_bkgds[0],bool(options.scale_bkgds[1]))

    groupsOff = options.turn_off_groups.split(',')
    if groupsOff[0] != '':
        for groupOff in groupsOff:
            print 'Turning group off the fit', groupOff
            sens.TurnGroupOff(groupOff)

    sens.GenAndFitData(options.number_runs,options.years,options.signal_counts,options.random_rate)
    # sens.~nEXOSensitivity()
    cmd = 'mv %s %s' % (outFileName,realPathDone)
    print(cmd)

    import time
    time.sleep(1)
    os.system(cmd)
    time.sleep(1)
    outfile = ROOT.TFile("%s/%s_%0.1f_years_%0.2f_counts_%04i.root"%(realPathDone,options.output_name,options.years,options.signal_counts,options.random_seed))
    outtree = outfile.Get("tree")

    ratiohist = ROOT.TH1F("ratiohist","ratiohist",100,0,10)
    outtree.Draw("nll_ratio>>ratiohist")
    xq=array.array('d',[.9])
    yq=array.array('d',[0])
    ratiohist.GetQuantiles(1,yq,xq)
    percent90=yq[0]
    print("magic number of %f counts: %f"%(options.signal_counts, percent90))