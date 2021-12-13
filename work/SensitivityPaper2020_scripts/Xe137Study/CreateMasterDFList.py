import pandas as pd
import numpy as np
import pickle
import os


critical_lambda_dir = '/p/lustre2/lenardo1/sensitivity_output/Mar7_Xe137_CriticalLambda_D024/'
scale_factor_tag = '_03.0x_'
num_iterations = 10

all_files = os.listdir( critical_lambda_dir )

h5_files = []


for i in range(100):
    #if i==15 or i==5: continue
    h5_files.append( [filename for filename in all_files if scale_factor_tag in filename 
                                                         and filename.endswith('_{}.h5'.format(i))] )

print('{} files found'.format(len(h5_files[0])))

dflist = []
inputcounts = []
sorteddflist = []

for i in range(num_iterations):
    
    print('Loading data for {}'.format(i))
    
    dflist.append([])
    inputcounts.append([])
    sorteddflist.append([])

    for filename in h5_files[i]:
        
        dflist[i].append( pd.read_hdf( critical_lambda_dir + filename ) )
        inputcounts[i].append( float(filename.split('_')[-3]) )
    sorteddflist[i] = [x for _,x in sorted( zip(inputcounts[i], dflist[i]) ) ]

print('Data loaded.') 

# First, combine the 0-9 lists to get one master list
masterdflist = []

for i in range(len(sorteddflist[0])):
    #masterdflist.append( sorteddflist1[i] )
    
    temp_list = [ sorteddflist[j][i] for j in range(num_iterations) ]
    masterdflist.append( pd.concat( temp_list, ignore_index=True ) )


with open(critical_lambda_dir + 'master_df_list_rn222{}D-024.pkl'.format(scale_factor_tag),'wb') as pklfile:
     pickle.dump(masterdflist,pklfile)

