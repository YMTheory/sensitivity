
import ROOT
import copy

ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')


filename = '../../v3/results/done/disc_hamamatsu_0nu_eff_tpc_rdm_5.0_years_0.0_counts_0001.root' #allfits_hamamatsu_0nu_eff_tpc_rdm_5.0_years_0.0_counts.root'

chain = ROOT.TChain('tree')
chain.Add(filename)

avgcorr = None
parnames = {}
tot = 0

max = 0

n = chain.GetEntries()
for i in range(n):
    chain.GetEntry(i)
    result = chain.fitres_sig
    pars = result.floatParsFinal()

    if not avgcorr:
        avgcorr = ROOT.TMatrixDSym(pars.getSize())
        for j in range(pars.getSize()):
            parnames[pars.at(j).GetName()] = j

    if chain.covQual_sig == 3 and chain.stat_sig == 0:
        tot += 1
        avgcorr += result.correlationMatrix()

        if abs(result.correlation('num_LXeBb0n','num_InternalU238')) > abs(max):
            max = result.correlation('num_LXeBb0n','num_InternalU238')
            print i, max

avgcorr *= 1./tot

avgcorr.Print()
        
print parnames    
    
print max
    
