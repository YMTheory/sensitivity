# Import useful libraries for analysis
import sys
import os
import pandas as pd
import histlite as hl
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

sys.path.append('../../modules')

# Import the nEXO sensitivity classes
import nEXOFitWorkspace
import nEXOFitModel
import nEXOFitLikelihood

# path_home = '/p/lustre2/czyz1/nexo_sensitivity/work'
# path_result = '/p/lustre2/nexouser/czyz1'
path_home = '/Users/czyz1/lc-home/nexo_sensitivity/work'
path_result = '/Users/czyz1/lc-nexouser'

base = "{}/czyz1/workdir/lambda/20_12_22_DNN1_023/".format(path_result)
# base = "/Users/czyz1/lc-nexouser/workdir/lambda/20_12_22_DNN1_023/"
resolution_list = ['0.008']#, '0.009', '0.01', '0.011', '0.012', '0.013', '0.014', '0.015', '0.016']
#resolution_list = ['0.013']
lambdas = {}
for resolution in resolution_list:
    lambdas[resolution] = []
    fixed_fit_covar = np.array([], dtype=bool)

    for run in range(122):
        file = base + 'critical_lambda_eres_{}_resolution_{}.h5'.format(run, resolution)
        if not os.path.exists(file):
            temp = 0
            for run_nest in [run-1, run+1]:
                file = base + 'critical_lambda_eres_{}_resolution_{}.h5'.format(run-1, resolution)
                data = pd.read_hdf(file, 'df')
                temp += np.percentile(data['lambda'], 90)/2

            lambdas[resolution].append(temp)
        
        else:

            data = pd.read_hdf(file, 'df')

            fixed_fit_covar = np.append(fixed_fit_covar, data['fixed_fit_covar'])

            lambdas[resolution].append(np.percentile(data['lambda'], 90))
        
    print(lambdas[resolution])


for resolution in resolution_list:
    with open(base + 'lambdas_res={}.txt'.format(resolution), 'w') as f:
        for item in lambdas[resolution]:
            f.write("%s\n" % item)
