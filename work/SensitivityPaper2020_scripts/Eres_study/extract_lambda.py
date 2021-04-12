# Import useful libraries for analysis
import sys
import os
import pandas as pd
import histlite as hl
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

sys.path.append('../../../modules')

# Import the nEXO sensitivity classes
import nEXOFitWorkspace
import nEXOFitModel
import nEXOFitLikelihood

path_home = '/p/lustre2/czyz1/nexo_sensitivity/work'
path_result = '/p/lustre2/nexouser/czyz1'
# path_home = '/Users/czyz1/lc-home/nexo_sensitivity/work'
# path_result = '/Users/czyz1/lc-nexouser'

date = '21_03_22'
db_list = ['024']

# resolution_list = ['0.008', '0.009', '0.01', '0.011', '0.012', '0.013', '0.014', '0.015', '0.016']
resolution_list = ['0.008']#, '0.011', '0.015', '0.018']
#resolution_list = ['0.013']
livetime_list = ['0.25']
lambdas = {}
for db in db_list:
    resolution = '0.008'
    # for resolution in resolution_list:
    for livetime in livetime_list:
        base = "{}/workdir/lambda/{}_DNN1_{}/".format(path_result, date, db)
        lambdas[resolution] = []
        fixed_fit_covar = np.array([], dtype=bool)

        for run in range(0, 122):
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

                # fixed_fit_covar = np.append(fixed_fit_covar, data['fixed_fit_covar'])

                # for ffc, lam in zip(data['fixed_fit_covar'], data['lambda']):
                #     if np.isnan(lam):
                #         ffc = False
                #
                # if np.percentile(np.array([a for a, b, c, d, e in zip(data['lambda'],
                #                                                    data['best_fit_converged'], data['best_fit_covar'],
                #                                                    data['fixed_fit_covar'],
                #                                                    data['fixed_fit_converged']) if
                #                         b and c and d and e and a >= 0]), 90) > 3.3:
                #     print(run)
                #     good = []
                #     too_big = []
                #     negative = []
                #     good_b = []
                #     too_big_b = []
                #     negative_b = []
                #     for index, lamb in enumerate(data['lambda']):
                #         if lamb >= 10:
                #             too_big.append([data['best_fit_parameters'].values[index]['Num_Xe137_and_Ar42']])
                #             too_big_b.append([data['best_fit_parameters'].values[index]['Num_B8nu']])
                #         elif lamb >= 0 and lamb < 10:
                #             good.append(data['best_fit_parameters'].values[index]['Num_Xe137_and_Ar42'])
                #             good_b.append([data['best_fit_parameters'].values[index]['Num_B8nu']])
                #         else:
                #             negative.append(data['best_fit_parameters'].values[index]['Num_Xe137_and_Ar42'])
                #             negative_b.append(data['best_fit_parameters'].values[index]['Num_B8nu'])
                #
                #     print('Mean too_big = {}, Mean too_big_b = {}, Median too_big = {}, Median too_big_b = {}'.format(
                #         np.mean(too_big), np.mean(too_big_b),  np.median(too_big), np.median(too_big_b)))
                #     print('Mean good = {}, Mean good_b = {}, Median good = {}, Median good_b = {}'.format(
                #         np.mean(good), np.mean(good_b), np.median(good), np.median(good_b)))
                #     print('Mean negative = {}, Mean negative_b = {}, Median negative = {}, Median negative_b = {}'.format(
                #         np.mean(negative), np.mean(negative_b), np.median(negative), np.median(negative_b)))
                #
                lambdas[resolution].append(np.percentile(np.array([a for a, b, c, d, e in zip(data['lambda'],
                                 data['best_fit_converged'], data['best_fit_covar'], data['fixed_fit_covar'],
                                 data['fixed_fit_converged']) if b and c and d and e and a>=-.01 and a<=10]), 90))

        print(lambdas[resolution])


    for resolution in resolution_list:
        with open(base + 'lambdas_res={}.txt'.format(resolution), 'w') as f:
            for item in lambdas[resolution]:
                f.write("%s\n" % item)


# i=-1
# color = ['b','g','r','k']
# for run in range(21,25):
#     i+= 1
#     file = base + 'critical_lambda_eres_{}_resolution_{}.h5'.format(run, resolution)
#     if not os.path.exists(file):
#         temp = 0
#         for run_nest in [run-1, run+1]:
#             file = base + 'critical_lambda_eres_{}_resolution_{}.h5'.format(run-1, resolution)
#             data = pd.read_hdf(file, 'df')
#             temp += np.percentile(data['lambda'], 90)/2
#         lambdas[resolution].append(temp)
#     else:
#         data = pd.read_hdf(file, 'df')
#
#     foo = np.array([a for a, b, c, d, e in zip(data['lambda'],
#                                                data['best_fit_converged'], data['best_fit_covar'],
#                                                data['fixed_fit_covar'], data['fixed_fit_converged']) if
#                     b and c and d and e])
#     n, bins, patches = plt.hist(foo, bins=len(range(-40, 40)), density=False, alpha=0.02, log=True)
#     plt.scatter(bins[:-1] + 0.5 * (bins[1:] - bins[:-1]), n, marker='o', c=color[i], s=10, alpha=1)