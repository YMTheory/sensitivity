
import ROOT

resultsFileName = 'results_5.0_years_0.0_counts.root'
plotsFileName = resultsFileName.replace('results','plots')

inFileName = '../results/full/%s' % (resultsFileName) 
outFileName = '../results/plots/%s' % (plotsFileName)

ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')

nu = ROOT.nEXOUtil()
nu.PlotBb0nUL(inFileName,outFileName)
