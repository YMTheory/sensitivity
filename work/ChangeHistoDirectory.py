
import ROOT
import sys


def main():

    if len(sys.argv) < 2:
        sys.exit('Usage: %s dir_name' % sys.argv[0])

    ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')

    # This is the individual_histos directory
    dirout = sys.argv[1]
    fin = ROOT.TFile.Open('../tables/Summary_v73_2016-09-09_0nu_allbkgs.root', 'read')
    tin = fin.Get('ExcelTableValues')

    table = ROOT.ExcelTableValues()

    tout = tin.CloneTree(0)
    tin.SetBranchAddress('table', table)
    tout.SetBranchAddress('table', table)

    print(tout.GetEntries())
    for i in range(tin.GetEntries()):
        tin.GetEntry(i)
        filename = table.fFileName
        dirin = filename[:filename.rfind('/')]
        print(table.fPdf, filename.replace(dirin, dirout))
        table.fFileName = filename.replace(dirin, dirout)

        tout.Fill()

    print(tout.GetEntries())

    fout = ROOT.TFile.Open('../tables/Summary_v73_2016-09-09_0nu_allbkgs_llnl.root', 'recreate')
    tout.Write()
    fout.Close()

    fin.Close()

if __name__ == "__main__":
    main()
