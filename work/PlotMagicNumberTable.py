
import ROOT
import os,sys
from optparse import OptionParser
import array
from os import listdir
from os.path import isfile, join
import numpy as np

#sample usage: python MakeMagicNumberTable.py -n 100 -r 1 -y 10.0 -c 0 -s 1 -d ../results -o recalctest_all -t ../tables/Summary_v73_2016-09-09_0nu_allbkgs_llnl.root -m 3 --ssfrac-improvement 1.0 --rn222-rate-correction 1.0 -M ../results/all_bkgs/


nknots = 9
spline_xn = array.array('d', [0.1, 1., 5., 7., 10., 15, 20, 30, 50.])
magic_numbers = []

def spline_func(x, par):

    xx = x[0]

    yn = array.array('d', [ par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8] ])

    sp3 = ROOT.TSpline3("sp3", spline_xn, yn, nknots, "b1e1", 0, 0)
    output = sp3.Eval(xx)
    del sp3

    return output


def main():

    usage = "usage: python RunSensitivity.py [options]"
    parser = OptionParser(usage)
    
    parser.add_option("-M","--magic-number-dir", nargs=1) # name of the output file
    options,args = parser.parse_args()

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

#    print(magic_numbers)
#    exit()
    # Magic Numbers need to be sorted by numcounts
    ndtemp = np.array(magic_numbers)
    np_magic = ndtemp[ndtemp[:, 0].argsort()]

    for c1, c2, c3 in np_magic:
        print c1, "\t", c2, "\t", c3

    gr = ROOT.TGraph(len(np_magic), array.array('d', np_magic[:, 0]), array.array('d', np_magic[:, 1]))
    gr.SetMarkerStyle(22)
    gr.SetMarkerColor(ROOT.kRed)
    #gr.Draw("AP")
    outfile = ROOT.TFile("critical_lambda_plot_file.root","recreate")
    gr.SetName("gr")
    gr.Write()
    outfile.Close()
    exit()

    # ROOT.gSystem.ProcessEvents()

    xmin = 0.
    xmax = 50.
    npars = 9

    f_spline4 = ROOT.TF1("f_spline", spline_func, xmin, xmax, npars)
    f_spline4.SetLineColor(ROOT.kBlue)
    f_spline4.SetParameters(3, 3, 3, 3, 3, 3, 3, 3, 3)
    # f_spline4.FixParameter(0, magic_numbers.front());
    # f_spline4.FixParameter(8, magic_numbers.back());

    gr.Fit(f_spline4, "R0")

    spline_yn = f_spline4.GetParameters()
    for i in range(0, nknots):
        print spline_xn[i], spline_yn[i]
        #sens.AddMagicNumber(numcounts,percent90)


if __name__ == "__main__":
    main()
