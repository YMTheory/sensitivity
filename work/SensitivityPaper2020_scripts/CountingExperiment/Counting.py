import numpy as np
import ROOT
import scipy.special as sp
import scipy.stats as st

Bb0nEfficiency = 0.5756

# see equation 5 in https://arxiv.org/pdf/2009.07249.pdf
def GetDiscoveryCCGV(Bkg, Signal, BkgError=0):
    if BkgError == 0: 
        return np.sqrt(2 * ((Signal + Bkg) * np.log((1 + Signal/Bkg)) - Signal ) )
    else: 
        E = Error * B
        FirstLog = ((S + B) * (B + E**2))/(B**2 + (S+B) * E**2)
        SecondLog = 1 + ((E**2 * S) / (B * (B + E**2)))
        Full = 2 * ((S + B) * np.log(FirstLog) - (B**2/E**2) * np.log(SecondLog) )
        return np.sqrt(Full)

# see equation 6 in https://arxiv.org/pdf/2009.07249.pdf
def GetSensitivityCCGV(Bkg, Signal, BkgError=0): 
    if BkgError == 0:
        return np.sqrt(2 * (Signal * np.log(1 + Signal/Bkg)) )
    else: 
        return 0

# see equation 8 in https://arxiv.org/pdf/2009.07249.pdf
# note that the expression is slightly different since scipy's implementaion
# of the incomplete lower gamma functions already includes the normalizaion
def GetDiscoveryAsimov(Bkg, Signal, BkgError=0):
    if BkgError == 0: 
        PValue = sp.gammainc(Signal + Bkg, Bkg)
    else: 
        PValue = sp.betainc(1)
    return np.abs(st.norm.ppf(PValue))

# see equation 9 in https://arxiv.org/pdf/2009.07249.pdf
# note that the expression is slightly different since scipy's implementaion
# of the incomplete upper gamma functions already includes the normalizaion
def GetSensitivityAsimov(Bkg, Signal, BkgError=0):
    PValue = sp.gammaincc(Bkg + 1, Signal + Bkg)
    return np.abs(st.norm.ppf(PValue))

def ConvertBkgIndexToBkgCounts(BkgIndex, Livetime, XenonMass):
    return BkgIndex * XenonMass * Livetime

# find the crossing, i.e. which signal number corresponds to a certain CL by 
# using a linear interpolation of significance vs signal count data
def GetCrossing(Signal, Significance, CL):
    XInterp = np.linspace(Signal[0], Signal[-1], 1000)
    Interpolation = np.interp(XInterp, Signal, Significance)
    Diff = np.abs(Interpolation - CL)
    Cut = np.where(Diff == np.min(Diff))[0][0]
    return XInterp[Cut]

def ConvertSignalToHalflife(Signal, Efficiency, XenonMass, Livetime):
    mmass_xe136 = 135.907214 # in g/mol 
    mmass_xe134 = 133.905393 # in g/mol 
    Avogadro = 6.022E23 # in 1/mol
    frac_136 = 0.9
    frac_134 = 0.1
    NXenonAtoms = XenonMass*1000.0*frac_136 / (frac_136*mmass_xe136 + frac_134*mmass_xe134) * Avogadro
    return Efficiency * NXenonAtoms * Livetime * np.log(2) / Signal

def RunCountingExperiment(XenonMass, Livetime, BkgIndex, Scaling):
    # Let's define the range over which we want to test different signal hypothesis. Here we choose a fine binning for average number of signal events 
    SignalRange = np.arange(0.01,1,0.01)
    SignalRange = np.append(SignalRange[:-1], np.arange(1,50,0.2))

    Significance = {}
    CLs = ['Discovery', 'Sensitivity']
    Methods = ['CCGV', 'Asimov']

    for cl in CLs: 
        Significance[cl] = {}
        for method in Methods: 
            Significance[cl][method] = {}
            Significance[cl][method]['Values'] = np.zeros((len(XenonMass), len(Livetime), len(SignalRange)))
            Significance[cl][method]['Crossing'] = np.zeros((len(XenonMass), len(Livetime)))

    # Scan over all fiducial masses that were included above
    for ii,Mass in enumerate(XenonMass): 

        # Scan over a fine range of livetimes for the experiment as defined above
        for jj,Time in enumerate(Livetime):

            # The number of background counts is given by the background index for the fiducial mass in question times the fiducial mass and the livetime of the experiment
            Bkg = BkgIndex[Mass] * Mass * Time / Scaling

            # Scan over a big range of signal hypothesis to find the number of signal counts that yield a certain confidence level
            for kk,Counts in enumerate(SignalRange):    

                # Calculate discovery significance based on different methods
                Significance['Discovery']['CCGV']['Values'][ii,jj,kk] = GetDiscoveryCCGV(Bkg=Bkg, Signal=Counts)
                Significance['Discovery']['Asimov']['Values'][ii,jj,kk] = GetDiscoveryAsimov(Bkg=Bkg, Signal=Counts)
                Significance['Sensitivity']['CCGV']['Values'][ii,jj,kk] = GetSensitivityCCGV(Bkg=Bkg, Signal=Counts)
                Significance['Sensitivity']['Asimov']['Values'][ii,jj,kk] = GetSensitivityAsimov(Bkg=Bkg, Signal=Counts)

                # Stop the signal hypothesis scan once all methods have reached a certain significance threshold
                if Significance['Discovery']['CCGV']['Values'][ii,jj,kk] > 3.2 and Significance['Discovery']['Asimov']['Values'][ii,jj,kk] > 3.2: 

                    # Find the number of signal events that correspond to the predefined CL by linear interpolation
                    Significance['Discovery']['CCGV']['Crossing'][ii,jj] = GetCrossing(SignalRange[:kk], Significance['Discovery']['CCGV']['Values'][ii,jj,:kk], CL=3.0)
                    Significance['Discovery']['Asimov']['Crossing'][ii,jj] = GetCrossing(SignalRange[:kk], Significance['Discovery']['Asimov']['Values'][ii,jj,:kk], CL=3.0)
                    Significance['Sensitivity']['CCGV']['Crossing'][ii,jj] = GetCrossing(SignalRange[:kk], Significance['Sensitivity']['CCGV']['Values'][ii,jj,:kk], CL=st.norm.ppf(0.90))
                    Significance['Sensitivity']['Asimov']['Crossing'][ii,jj] = GetCrossing(SignalRange[:kk], Significance['Sensitivity']['Asimov']['Values'][ii,jj,:kk], CL=st.norm.ppf(0.90))

                    X,Y = np.meshgrid(Livetime, XenonMass)

                    Significance['Discovery']['CCGV']['Half-life'] = ConvertSignalToHalflife(Signal=Significance['Discovery']['CCGV']['Crossing'], Efficiency=Bb0nEfficiency, XenonMass=Y, Livetime=X) 
                    Significance['Discovery']['Asimov']['Half-life'] = ConvertSignalToHalflife(Signal=Significance['Discovery']['Asimov']['Crossing'], Efficiency=Bb0nEfficiency, XenonMass=Y, Livetime=X) 
                    Significance['Sensitivity']['CCGV']['Half-life'] = ConvertSignalToHalflife(Signal=Significance['Sensitivity']['CCGV']['Crossing'], Efficiency=Bb0nEfficiency, XenonMass=Y, Livetime=X) 
                    Significance['Sensitivity']['Asimov']['Half-life'] = ConvertSignalToHalflife(Signal=Significance['Sensitivity']['Asimov']['Crossing'], Efficiency=Bb0nEfficiency, XenonMass=Y, Livetime=X) 
                    break
    
    return Significance