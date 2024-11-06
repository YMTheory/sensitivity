from os import WSTOPSIG
import numpy as np

from neutrino_source import neutrino_source
from detector import detector
from oscillation import *
from xsection import xsection

## This script will calculate the weight factor for experiments normalized to the nonimal nEXO configuration.
## One way is to save them all into a json files 
## Another way is to calculate it real time

def weight(source, det, int_type='ES'):

    # The geometry effect (like different detector geometrys and distances will be considered in the Monte Carlo sampling) .
    Enu = source.energies[0]
    epsilon = det.efficiency
    enrichment = det.enrichment
    # We will only consider the following terms:
    # 1. cross section
    if int_type == 'CC':
        w_xsec = 1.0
    elif int_type == 'ES':
        xsec = xsection()
        w_xsec = xsec.total_xsec_ES(Enu) / xsec.total_xsec_CC(Enu) 
    # 2. enrichment:
    w_enrichment = enrichment / 0.90
    # 3. epsilon:
    w_efficiency = epsilon / 1.0
    # 4. total flux ratio:
    flux = source.flux_integral()
    flux0 = source.flux_integral_nominal()
    w_flux = flux / flux0

    return w_xsec * w_enrichment * w_efficiency * w_flux


def flux_activity(activity0, t, halflife=27.7):
    return activity0 * np.exp(-t*(0.693/halflife))
    
def flux_integral(activity0, t, halflife=27.7):
    lbd = 0.693 / halflife
    return -1/lbd*(flux_activity(activity0, t, halflife=halflife) - flux_activity(activity0, 0., halflife=halflife))        
        

def weight_exps(int_types, efficiencies, enrichments, activities, times, volumes, Enu=0.75):
    # We need to give a list for each of these arguments, where the 1st is the detector we want to calculate, the 2nd is the one we use to normalize.
    xsecs = []
    xs = xsection()
    for ty in int_types:
        if ty == 'ES':
            xsecs.append( xs.total_xsec_ES(Enu) )
        elif ty == 'CC':
            xsecs.append( xs.total_xsec_CC(Enu) )
    w_xsec = xsecs[0] / xsecs[1]

    w_efficiency = efficiencies[0] / efficiencies[1]
    
    w_enrichment = enrichments[0] / enrichments[1]

    flux0 = flux_integral(activities[0], times[0])
    flux1 = flux_integral(activities[1], times[1])
    w_flux = flux0 / flux1

    w_volume = volumes[0] / volumes[1]

    return w_xsec * w_efficiency * w_enrichment * w_flux * w_volume
     


