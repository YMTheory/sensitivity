# Import useful libraries for analysis
import sys
import os
import pandas as pd
import numpy as np

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
db = '024'

resolution = '0.008'
lambdas = {}
livetime_list = ['0.5', '1', '2', '5']

for livetime in livetime_list:

    base = "{}/workdir/lambda/{}_DNN1_{}/".format(path_result, date, db)
    lambdas = [] #np.array([])
    fixed_fit_covar = np.array([], dtype=bool)

    for iteration in range(122):

        sub_lam = np.array([])
        for subit in range(5):

            file = base + 'critical_lambda_iteration_{}_lt_years_{}_numits={}.h5'.format(iteration, livetime, subit)
            print('Livetime = {}, iteration = {}, subit = {}'.format(livetime, iteration, subit))
            if not os.path.exists(file):
                pass
                # temp = 0
                # for run_nest in [it-1, it+1]:
                #     file = base + 'critical_lambda_eres_0_resolution_{}_numits={}.h5'.format(resolution, run_nest)
                #     data = pd.read_hdf(file, 'df')
                #     temp += np.array([a for a, b, c, d, e in zip(data['lambda'],
                #                  data['best_fit_converged'], data['best_fit_covar'], data['fixed_fit_covar'],
                #                  data['fixed_fit_converged']) if b and c and d and e and a>=0 and a<=10])
                #
                # lambdas.append(temp)

            else:
                data = pd.read_hdf(file, 'df')
                sub_lam = np.append(sub_lam, np.array([a for a, b, c, d, e in zip(data['lambda'],
                                             data['best_fit_converged'], data['best_fit_covar'],
                                             data['fixed_fit_covar'],
                                             data['fixed_fit_converged']) if b and c and d and e and a >= -.001 and a <= 10]))
                # lambdas = np.append(lambdas, np.array([a for a, b, c, d, e in zip(data['lambda'],
                #                              data['best_fit_converged'], data['best_fit_covar'],
                #                              data['fixed_fit_covar'],
                #                              data['fixed_fit_converged']) if b and c and d and e and a >= -.001 and a <= 10]))

        lambdas.append(np.percentile(sub_lam, 90))

    with open(base + 'lambdas_lt={}.txt'.format(livetime), 'w') as f:
        for item in lambdas: #[resolution]:
            f.write("%s\n" % item)

        # print(lambdas)
        # print(np.percentile(lambdas, 90))


# with open(base + 'lambdas_res={}.txt'.format(resolution), 'w') as f:
#     for item in lambdas[resolution]:
#         f.write("%s\n" % item)


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