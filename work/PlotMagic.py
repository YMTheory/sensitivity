import ROOT
from os import listdir
from os.path import isfile, join
import sys
import array
import numpy as np

magic_number_dir = "../results/all_bkgs/"

magic_numbers = []


def main():
    magic_files = {f[:-9] for f in listdir(magic_number_dir) if isfile(join(magic_number_dir, f))}
    # print(magic_files)

    # c1 = ROOT.TCanvas( 'c1', 'c1')
    for f in magic_files:
        numcounts = float(f.split('_')[-3])
        chain = ROOT.TChain("tree")
        chain.Add(join(magic_number_dir, f) + "*")
        print("Processing ", f)
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
        print(c1, "\t", c2, "\t", c3)

    inter = ROOT.Math.Interpolator(10, ROOT.Math.Interpolation.kCSPLINE)
    inter.SetData(len(np_magic), array.array('d', np_magic[:, 0]), array.array('d', np_magic[:, 1]))

    x = np.arange(0, max(np_magic[:, 0]), 0.1)
    y = np.array(map(inter.Eval, x))

    c2 = ROOT.TCanvas('c2', 'c2')
    gr = ROOT.TGraph(len(np_magic), array.array('d', np_magic[:, 0]), array.array('d', np_magic[:, 1]))
    gr.SetMarkerStyle(22)
    gr.Draw("AP")

    gr2 = ROOT.TGraph(len(y), x, y)
    gr2.Draw("SAME L")
    ROOT.gSystem.ProcessEvents()
    c2.Modified()
    c2.Update()
    print("press a key:")
    sys.stdin.read(1)

if __name__ == "__main__":
    main()
