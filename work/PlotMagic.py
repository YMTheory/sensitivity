import ROOT
from os import listdir
from os.path import isfile, join
import sys
import array
import numpy as np
import itertools

# FIXME: rewrite as a class of tools to fit and plot lamdba critical

magic_number_dir = "../results/all_bkgs-1DoF/"
scan_dir = "../results-0signal/done/"

# magic_number_dir = "../results-BaTag/done/"
# scan_dir = "../results-BaTag-0signal/done/"

plot_scan = False
plot_scan2 = True

# position of the knots that form the spline used to approximate the lambda critical curve from the data
# these might need to be tweaked depending on the actual shape of the curve
spline_xn = array.array('d', [0.1, 1., 5., 7., 10., 15, 20, 30, 50.])  # all_bkgs
# spline_xn = array.array('d', [0.1, 0.5, 1., 2., 3.5, 4.5, 10, 20, 50.])   # Ba_tag
nknots = len(spline_xn)

magic_numbers = []


def spline_func(x, par):
    xx = x[0]

    yn = array.array('d', [par[i] for i in range(0, nknots)])

    sp3 = ROOT.TSpline3("sp3", spline_xn, yn, nknots, "b1e1", 0, 0)
    output = sp3.Eval(xx)
    del sp3

    return output


def plot_magic():
    h2 = ROOT.TH2F("h2", "h2;Signal Hypothesis;Critical Lambda", 100, 0., spline_xn[-1], 100, 0, 5)
    if plot_scan:
        scan_files = {f[:-9] for f in listdir(scan_dir) if isfile(join(scan_dir, f))}

        for f in scan_files:
            numcounts = float(f.split('_')[-3])
            chain = ROOT.TChain("tree")
            chain.Add(join(scan_dir, f) + "*")
            print "Processing ", join(scan_dir, f)
            chain.Draw("nll_ratio:" + str(numcounts) + ">>+h2", "nll_ratio>=0", "goff")

    magic_files = {f[:-9] for f in listdir(magic_number_dir) if isfile(join(magic_number_dir, f))}
    # print(magic_files)

    # c1 = ROOT.TCanvas( 'c1', 'c1')
    for f in magic_files:
        numcounts = float(f.split('_')[-3])
        chain = ROOT.TChain("tree")
        chain.Add(join(magic_number_dir, f) + "*")
        print "Processing ", join(magic_number_dir, f)
        ratiohist = ROOT.TH1F("ratiohist", f, 400, 0, 20)
        chain.Draw("nll_ratio>>ratiohist",
                   "nll_ratio>=0 && stat_sig==0 && stat_bkg==0 && covQual_sig==3 && covQual_bkg==3", "goff")
        xq = array.array('d', [.9])
        yq = array.array('d', [0])
        ratiohist.GetQuantiles(1, yq, xq)
        percent90 = yq[0]
        # print numcounts, "\t", percent90, "\t",
        magic_numbers.append((numcounts, percent90 / 2., ratiohist.GetEntries()))

        # refresh the canvases
        # ROOT.gSystem.ProcessEvents()
        # c1.Modified()
        # c1.Update()
        # print "press a key:"
        # sys.stdin.read(1)

        del ratiohist

    ndtemp = np.array(magic_numbers)
    np_magic = ndtemp[ndtemp[:, 0].argsort()]

    for c1, c2, c3 in np_magic:
        print c1, "\t", c2, "\t", c3

    l = ROOT.TLine(0, 1.35, spline_xn[-1], 1.35)
    l.SetLineColor(ROOT.kRed)

    mg = ROOT.TMultiGraph()

    c2 = ROOT.TCanvas('c2', 'c2')
    gr = ROOT.TGraph(len(np_magic), array.array('d', np_magic[:, 0]), array.array('d', np_magic[:, 1]))
    gr.SetMarkerStyle(22)
    gr.SetMarkerColor(ROOT.kBlue)
    gr.SetLineColor(ROOT.kBlue)
    # gr.Draw("AP")

    # gr2 = ROOT.TGraph(len(y), x, y)
    # gr2.Draw("SAME L")
    ROOT.gSystem.ProcessEvents()

    f_spline4 = ROOT.TF1("f_spline", spline_func, 0., spline_xn[-1], nknots)
    f_spline4.SetLineColor(ROOT.kBlue)
    map(f_spline4.SetParameter, range(0, nknots), itertools.repeat(2, nknots))

    gr.Fit(f_spline4, "R")

    spline_yn = f_spline4.GetParameters()
    for i in range(0, nknots):
        print spline_xn[i], spline_yn[i]

    if (plot_scan2):
        scan_files = {f[:-9] for f in listdir(scan_dir) if isfile(join(scan_dir, f)) and f[-9:-5] == '0001'}

        # nll_ratio = ROOT.TLeaf()
        for i in range(0, 50):
            nll_ratios = []
            print "Entry ", i
            for f in scan_files:
                numcounts = float(f.split('_')[-3])
                chain = ROOT.TChain("tree")
                chain.Add(join(scan_dir, f) + "0001*")
                nll_ratio = chain.GetLeaf("nll_ratio")
                covQual_sig = chain.GetLeaf("covQual_sig")
                chain.GetEntry(i)
                if covQual_sig.GetValue() != 3:
                    nll_ratios = []
                    break
                nll_ratios.append((numcounts, nll_ratio.GetValue() / 2.))

            if not nll_ratios: continue
            ndtemp = np.array(nll_ratios)
            np_nll_ratios = ndtemp[ndtemp[:, 0].argsort()]
            gr_tmp = ROOT.TGraph(len(np_nll_ratios), array.array('d', np_nll_ratios[:, 0]),
                                 array.array('d', np_nll_ratios[:, 1]))
            # gr_tmp.Fit("pol4", "q")
            mg.Add(gr_tmp)
            # gr_tmp.Draw("SAME L")

    mg.Add(gr)
    mg.Draw("A L")
    l.Draw()
    ROOT.gSystem.ProcessEvents()

    c3 = ROOT.TCanvas('c3', 'c3')
    h2.Draw("COLZ")
    gr.Draw("SAME LP")
    l.Draw()

    c2.cd()
    c2.Modified()
    c2.Update()
    print("press a key:")
    sys.stdin.read(1)


if __name__ == "__main__":
    plot_magic()
