import ROOT
import os
from optparse import OptionParser
import array
from os import listdir
from os.path import isfile, join

import itertools
import numpy as np

# sample usage: python MakeMagicNumberTable.py -n 100 -r 1 -y 10.0 -c 0 -s 1 -d ../results -o recalctest_all -t ../tables/Summary_v73_2016-09-09_0nu_allbkgs_llnl.root -m 3 --ssfrac-improvement 1.0 --rn222-rate-correction 1.0 -M ../results/all_bkgs/

# FIXME: several parts of this code should be imported from PlotMagic.py instead than replicating code here

# position of the knots that form the spline used to approximate the lambda critical curve from the data
# these might need to be tweaked depending on the actual shape of the curve
# use PlotMagic.py to find optimal values
spline_xn = array.array('d', [0.1, 1., 5., 7., 10., 15, 20, 30, 50.])
nknots = len(spline_xn)

magic_numbers = []


def spline_func(x, par):
    xx = x[0]

    yn = array.array('d', [par[i] for i in range(0, nknots)])

    sp3 = ROOT.TSpline3("sp3", spline_xn, yn, nknots, "b1e1", 0, 0)
    output = sp3.Eval(xx)
    del sp3

    return output


def main():
    libRelPath = '../lib/libnEXOSensitivity.so'
    ROOT.gSystem.Load(libRelPath)
    # ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

    usage = "usage: python RunSensitivity.py [options]"
    parser = OptionParser(usage)

    parser.add_option("--binning", nargs=6,
                      type=float)  # choose energy and SD binning, format: n-energy-bins, energy-min, energy-max, n-sd-bins, sd-min, sd-max

    parser.add_option("--ba-tag", action="store_true", default=False)  # run in Ba-tag scenario (only bb2n background)
    parser.add_option("--scale-bkgds", nargs=2,
                      type=float)  # scale all background counts, 2nd argument whether scale bb2n or not

    parser.add_option("-n", "--number-runs", nargs=1, type=int, default=1)  # number of toy fits to be produced
    parser.add_option("-s", "--random-seed", nargs=1, type=int, default=1)  # seed for randomizations
    parser.add_option("-y", "--years", nargs=1, type=float, default=5.)  # detector livetime in years
    parser.add_option("-t", "--tree",
                      nargs=1)  # ROOT file with tree for input counts expectation (produced from summary table)
    parser.add_option("-d", "--output-dir", nargs=1)  # name of the output directory
    parser.add_option("-o", "--output-name", nargs=1)  # name of the output file
    parser.add_option("-M", "--magic-number-dir", nargs=1)  # name of the output file

    parser.add_option("-c", "--signal-counts", nargs=1, type=float,
                      default=0.0)  # included signal counts (for discovery potential)

    parser.add_option("-r", "--random-rate", nargs=1, type=int,
                      default=0)  # rate (in evt^-1) at which the activity is randomized
    parser.add_option("-m", "--method-bkgd", nargs=1, type=int,
                      default=ROOT.nEXOSensitivity.kRdmCV)  # method to choose activities : ROOT.nEXOSensitivity.kRdmCV, kUL, kPosUL

    parser.add_option("-g", "--group", nargs=1)  # modify counts in group given by name here
    parser.add_option("-v", "--bkgd-counts", nargs=1,
                      type=float)  # modify to counts given here for group given by name above
    parser.add_option("--ssfrac-improvement", nargs=1, type=float,
                      default=1.)  # change of SS/MS fraction for gamma-like PDFs
    parser.add_option("--rn222-rate-correction", nargs=1, type=float, default=1.)  # change of Rn222 rate

    parser.add_option("--turn-off-groups", nargs=1, type=str,
                      default='')  # turn these groups off the fit (set zero expectation counts and do not include in fit)

    options, args = parser.parse_args()
    print 'Using options:', options

    realPath = os.path.realpath(options.output_dir)
    realPathWork = realPath + '/working'
    realPathDone = realPath + '/done'
    os.system('mkdir -p %s' % (realPathWork))
    os.system('mkdir -p %s' % (realPathDone))
    outFileName = '%s/%s_%0.1f_years_%0.1f_counts_%04i.root' % (
        realPathWork, options.output_name, options.years, options.signal_counts, options.random_seed)

    sens = ROOT.nEXOSensitivity(options.random_seed, options.tree)
    sens.fResultFileName = outFileName

    sens.fRunBkgdOnlyFit = False
    sens.fExpectCountMethod = options.method_bkgd  # ROOT.nEXOSensitivity.kRdmCV #kUL #kPosUL #kRdmCV

    # The critical lambda is extracted from all files in the magic_number_dir that have different random seeds
    magic_files = {f[:-9] for f in listdir(options.magic_number_dir) if isfile(join(options.magic_number_dir, f))}
    # print(magic_files)

    # c1 = ROOT.TCanvas( 'c1', 'c1')
    for f in magic_files:
        numcounts = float(f.split('_')[-3])
        chain = ROOT.TChain("tree")
        chain.Add(join(options.magic_number_dir, f) + "*")
        print "Processing ", join(options.magic_number_dir, f)
        ratiohist = ROOT.TH1F("ratiohist", f, 400, 0, 20)
        chain.Draw("nll_ratio>>ratiohist",
                   "nll_ratio>=0 && stat_sig==0 && stat_bkg==0 && covQual_sig==3 && covQual_bkg==3", "goff")
        xq = array.array('d', [.9])
        yq = array.array('d', [0])
        ratiohist.GetQuantiles(1, yq, xq)
        percent90 = yq[0]
        # print numcounts, "\t", percent90, "\t",
        magic_numbers.append((numcounts, percent90, ratiohist.GetEntries()))

        # refresh the canvases
        # ROOT.gSystem.ProcessEvents()
        # c1.Modified()
        # c1.Update()
        # print "press a key:"
        # sys.stdin.read(1)

        del ratiohist

    # Magic Numbers need to be sorted by numcounts
    ndtemp = np.array(magic_numbers)
    np_magic = ndtemp[ndtemp[:, 0].argsort()]

    for c1, c2, c3 in np_magic:
        print c1, "\t", c2, "\t", c3

    gr = ROOT.TGraph(len(np_magic), array.array('d', np_magic[:, 0]), array.array('d', np_magic[:, 1]))
    gr.SetMarkerStyle(22)
    gr.SetMarkerColor(ROOT.kRed)
    # gr.Draw("AP")
    # ROOT.gSystem.ProcessEvents()

    # Create the spline to be fitted to the data
    f_spline4 = ROOT.TF1("f_spline", spline_func, 0., spline_xn[-1], nknots)
    f_spline4.SetLineColor(ROOT.kBlue)
    map(f_spline4.SetParameter, range(0, nknots), itertools.repeat(2, nknots))
    # do the actual fit
    gr.Fit(f_spline4, "R0")

    # pass the parameters to the sensitivity object
    spline_yn = f_spline4.GetParameters()
    for i in range(0, nknots):
        print spline_xn[i], spline_yn[i]
        # sens.AddMagicNumber(numcounts,percent90)
        sens.AddMagicNumber(spline_xn[i], spline_yn[i])

    if options.group:
        sens.AddUserMeanCounts(options.group, options.bkgd_counts)

    sens.fVerboseLevel = 0
    sens.fSSFracImprovement = options.ssfrac_improvement
    sens.fRn222RateCorrection = options.rn222_rate_correction

    if options.binning:
        sens.SetBinning(int(options.binning[0]), float(options.binning[1]), float(options.binning[2]),
                        int(options.binning[3]), float(options.binning[4]), float(options.binning[5]))

    sens.SetBaTag(options.ba_tag)
    if options.scale_bkgds:
        sens.ScaleMeanBackgrounds(options.scale_bkgds[0], bool(options.scale_bkgds[1]))

    groupsOff = options.turn_off_groups.split(',')
    if groupsOff[0] != '':
        for groupOff in groupsOff:
            print 'Turning group off the fit', groupOff
            sens.TurnGroupOff(groupOff)

    sens.GenAndFitData(options.number_runs, options.years, options.signal_counts, options.random_rate)

    cmd = 'mv %s %s' % (outFileName, realPathDone)
    print(cmd)
    os.system(cmd)


if __name__ == "__main__":
    main()
