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

for file in os.listdir(pathToROOTPDFs):
    num_processed+=1
    if '.root' not in file:
        continue
    print('Loading {} at {:.4} seconds...\t({}/{})'.format(file,\
                                                        time.time()-start_time,\
                                                        num_processed,\
                                                        num_rootfiles))
    thisfile = uproot.open( (pathToROOTPDFs + file) )
    h_StandoffVsEnergySS_Smear = thisfile['h_StandoffVsEnergySS_Smear;1'].numpy()
    h_StandoffVsEnergyMS_Smear = thisfile['h_StandoffVsEnergyMS_Smear;1'].numpy()
    hh = hl.Hist( [ [0,1,2],\
                    h_StandoffVsEnergySS_Smear[1][0][0],\
                    h_StandoffVsEnergySS_Smear[1][0][1] ],\
                  [ h_StandoffVsEnergySS_Smear[0],\
                    h_StandoffVsEnergyMS_Smear[0] ] )
    hh = hh.rebin(1,h_StandoffVsEnergySS_Smear[1][0][0][0::10])
    hh = hh.rebin(2,h_StandoffVsEnergySS_Smear[1][0][1][0::26])
    
    thisrow = {'Filename':file, 'PDF':hh}
    rowslist.append(thisrow)
    
df_pdf = pd.DataFrame(rowslist)
df_pdf.to_hdf(outputFile,key='SimulationPDFs')

end_time = time.time()

print('Elapsed time = {} seconds ({} minutes).'.format( end_time-start_time, (end_time-start_time)/60. ) )
