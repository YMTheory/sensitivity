
import ROOT
import array

libRelPath = '../lib/libnEXOSensitivity.so'
ROOT.gSystem.Load(libRelPath)

filenames = '../results/full/alldisc_hamamatsu_0nu_eff_tpc_rdm_5.0_years_%0.1f_counts.root'
#filenames = '/data/data033/exo/software/nEXO_Sensitivity/from_marius/rep_Caio_discnsac_Ba/results/done/results_5.0_yrs_%0.1f_counts_*.root'
counts = array.array('d',range(11))

disc = ROOT.nEXOUtils.GetDiscoveryCounts(0.9,3, len(counts),counts, filenames)
print disc
