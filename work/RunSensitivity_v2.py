
import ROOT
import os,sys
from optparse import OptionParser

if __name__ == "__main__":

    usage = "usage: python RunFitting.py -m module -a action -c commands [-s]"
    parser = OptionParser(usage)
    
    parser.add_option("-n","--number-runs", nargs=1,type=int,default=1)
    parser.add_option("-m","--random-rate", nargs=1,type=int,default=0)
    parser.add_option("-r","--random-seed", nargs=1,type=int,default=1)
    parser.add_option("-y","--years", nargs=1,type=float,default=5.)
    parser.add_option("-c","--signal-counts", nargs=1,type=float,default=0.0)
    parser.add_option("-d","--output-dir", nargs=1)
    parser.add_option("-o","--output-name", nargs=1)
    parser.add_option("-t","--tree", nargs=1)

    parser.add_option("-g","--group", nargs=1)
    parser.add_option("-v","--bkgd-counts", nargs=1,type=float)
    

    options,args = parser.parse_args()   
    print 'Using options:', options

    libRelPath = '../lib/libnEXOSensitivity.so'
    ROOT.gSystem.Load(libRelPath)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

    realPath = os.path.realpath(options.output_dir)
    realPathWork = realPath + '/working'
    realPathDone = realPath + '/done'
    os.system('mkdir -p %s' % (realPathWork))
    os.system('mkdir -p %s' % (realPathDone))
    outFileName = '%s/%s_%0.1f_years_%0.1f_counts_%04i.root' % (realPathWork,options.output_name,options.years,options.signal_counts,options.random_seed)

    sens = ROOT.nEXOSensitivity(options.random_seed,options.tree)
    sens.fResultFileName = outFileName
    sens.fExpectCountMethod = ROOT.nEXOSensitivity.kUL #kRdmCV

    if options.group:
        sens.AddUserMeanCounts(options.group,options.bkgd_counts)
    sens.GenAndFitData(options.number_runs,options.years,options.signal_counts,options.random_rate)

    cmd = 'mv %s %s' % (outFileName,realPathDone)
    print cmd
    os.system(cmd)
