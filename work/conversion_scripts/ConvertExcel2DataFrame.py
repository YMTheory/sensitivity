########################################################################
##
## Main script to convert values in summary Excel table into ROOT tree
## Usage: 
## $python ConvertExcel2Root.py
##
########################################################################


import sys
import pandas as pd
import histlite as hl
import numpy as np
import NameDict
import time
import os
#import openpyxl
import nEXOExcelTableReader

#######################################################################################
# Input to the script

if len(sys.argv) == 4:
  inTableName = sys.argv[1] # '../tables/Summary_v68_2016-06-21_0nu.xlsx'
  outTableName = sys.argv[2] # '../tables/Summary_v68_2016-06-21_0nu.h5'
  pathToPDFs = sys.argv[3]
  if not os.path.exists(pathToPDFs):
     sys.exit('\nERROR: path to PDF files does not exist\n')
else:
  print('\n\nERROR: ConvertExcel2DataFrame.py requires 3 arguments');
  print('Usage:')
  print('\tpython ConvertExcel2DataFrame.py </path/to/inputExcelTable> <outputHDF5FileName> </path/to/hdf5PdfsFile>')
  sys.exit('\n')
   

#######################################################################################
# Script actually starts here

start_time = time.time()

excelTable = nEXOExcelTableReader.nEXOExcelTableReader(inTableName,pathToPDFs,config='./config/TUTORIAL_config.yaml') #excelTable.Print()

try: 
   excelTable.ConvertExcel2DataFrame()
except KeyError:
   sys.exit()
   
print( excelTable.components )
excelTable.components.to_hdf(outTableName,key='Components')

end_time = time.time()
print('Elapsed time = {} seconds ({} minutes).'.format( end_time-start_time, (end_time-start_time)/60. ) )
