import pandas as pd
import uproot
import histlite as hl
import numpy as np
import os
import sys
import time


if len(sys.argv) == 3:
    pathToROOTPDFs = sys.argv[1]
    outputFile = sys.argv[2]
    if not os.path.exists(pathToROOTPDFs):
           sys.exit('\nERROR: path to PDF files does not exist\n')
else:
        print('\n\nERROR: ConvertRootPDFsToHistlite.py requires 2 arguments');
        print('Usage:')
        print('\tpython ConvertRootPDFsToHistlite.py </directory/with/rootPDFs/> <outputHDF5FileName>')
        sys.exit('\n')

start_time = time.time()

num_rootfiles = len(os.listdir(pathToROOTPDFs))

rowslist = []
num_processed = 0

for filename in os.listdir(pathToROOTPDFs):
    num_processed+=1
    if '.root' not in filename:
        continue
    print('Loading {} at {:.4} seconds...\t({}/{})'.format(filename,\
                                                        time.time()-start_time,\
                                                        num_processed,\
                                                        num_rootfiles))
    thisfile = uproot.open( (pathToROOTPDFs + '/' + filename) )
    h_StandoffVsEnergySS_Smear = thisfile['h_StandoffVsEnergySS_Smear;1'].numpy()
    h_StandoffVsEnergyMS_Smear = thisfile['h_StandoffVsEnergyMS_Smear;1'].numpy()
    #print(len(h_StandoffVsEnergySS_Smear[1]))
    #print(len(h_StandoffVsEnergySS_Smear[0]))
    hh = hl.Hist( [ [0,1,2],\
                    h_StandoffVsEnergySS_Smear[1][0][0][60:],\
                    h_StandoffVsEnergySS_Smear[1][0][1][1:] ],\
                  [ h_StandoffVsEnergySS_Smear[0][60:,1:].astype(float),\
                    h_StandoffVsEnergyMS_Smear[0][60:,1:].astype(float) ] )
    #print(hh)
    #print(len(h_StandoffVsEnergySS_Smear[1][0][0][60::2]))
    hh = hh.rebin(1,h_StandoffVsEnergySS_Smear[1][0][0][60::2])
    hh = hh.rebin(2,h_StandoffVsEnergySS_Smear[1][0][1][1::2])
    
    thisrow = {'Filename':filename, 'Histogram':hh, 'HistogramAxisNames':['SS/MS','Energy (keV)','Standoff (mm)']}
    rowslist.append(thisrow)
    
df_pdf = pd.DataFrame(rowslist)
df_pdf.to_hdf(outputFile,key='SimulationHistograms')

end_time = time.time()

print('Elapsed time = {} seconds ({} minutes).'.format( end_time-start_time, (end_time-start_time)/60. ) )
