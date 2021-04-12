# Import sys, then tell python where to find the nEXO-specific classes
import sys
import os

sys.path.append('../../../modules')

######################################################################
# Load the arguments:
# if len(sys.argv) != 4:
#     print('\nERROR: incorrect number of arguments.\n')
#     print('Usage:')
#     print('\tpython Compute90PercentLimit_PythonCode.py ' + \
#           '<input_data_file> <critical_lamda_file> ' + \
#           '<output_dir>\n\n')
#     sys.exit()
#
# input_data_file = sys.argv[1]
# critical_lambda_file = sys.argv[2]
# output_dir = sys.argv[3]
#

#####################################################################

##########################################################################
##########################################################################
def FindIntersectionByQuadraticInterpolation(xvals, yvals, SplineFunction):
    if len(xvals) != len(yvals):
        print('ERROR: need same length arrays for ' + \
              'quadratic interpolation')
        raise ValueError
    # print('Xvals:')
    # print(xvals)
    # print('Yvals:')
    # print(yvals)

    # First, select only lambda values on the upward slope, so we find
    # the upper limit
    mask = np.zeros(len(yvals), dtype=bool)
    mask[1:] = (yvals[1:] - yvals[:-1]) > 0.
    # Next, select only values near the critical lambda threshold (~2.7)
    mask = mask & (yvals > 0.5) & (yvals < 6.)

    xfit = np.linspace(0., 100., 10000)

    if len(xvals[mask]) > 0:
        try:
            p = np.polyfit(xvals[mask], yvals[mask], 2)
            # print('Quadratic fit: {}'.format(p))
            yfit = p[0] * xfit ** 2 + p[1] * xfit + p[2]
        except np.RankWarning:
            p = np.polyfit(xvals[mask], yvals[mask], 1)
            # print('Linear fit: {}'.format(p))
            yfit = p[0] * xfit + p[1]
        yspline = SplineFunction(xfit)
        crossing_idx = np.where((yfit - yspline) > 0.)[0][0]
        crossing = xfit[crossing_idx]
    else:
        yfit = np.zeros(len(xfit))
        crossing_idx = -1
        crossing = -1.

    return xfit, yfit, crossing, crossing_idx


##########################################################################
##########################################################################
#
# # input_data_file = sys.argv[1]
# # critical_lambda_file = sys.argv[2]
# # output_dir = sys.argv[3]
#
# # Import useful libraries for analysis
# import pandas as pd
# import histlite as hl
# import numpy as np
# from matplotlib import pyplot as plt
# from scipy.interpolate import LSQUnivariateSpline
# from scipy.interpolate import UnivariateSpline
#
# # Import the nEXO sensitivity classes
# import nEXOFitWorkspace
# import nEXOFitModel
# import nEXOFitLikelihood
#
# dnn = 'DNN1'
# materialdb = "024"
# dates = ["21_02_08_{}_{}".format(dnn, materialdb)]#, "21_01_21_{}_023".format(dnn)]
# materialdb = "025"
# lamdate = "21_01_21_{}_{}".format(dnn, materialdb)
# # resolutions = ['0.008', '0.009', '0.01', '0.011', '0.012', '0.013', '0.014', '0.015', '0.016', '0.017', '0.018']
# # resolutions = ['0.011','0.015','0.018']
# resolutions = ['0.008']
# datasets = 100
#
# # path_home = '/Users/czyz1/lc-home/nexo_sensitivity/work'
# # path_result = '/p/lustre2/nexouser/czyz1'
# path_result = '/Users/czyz1/lc-nexouser'
#
# root_dir = '{}/output/h5/'.format(path_result)
# # lam_file_reses = ['0.011', '0.015', '0.018']
# lam_file_reses = ['0.008']
#
# # input_data_file = sys.argv[1]
# # critical_lambda_file = sys.argv[2]
# # output_dir = sys.argv[3]
#
# for index, resolution in enumerate(resolutions):
#
#     if resolution == '0.008':
#         iterlams = ['0.008']
#     else:
#         iterlams = ['0.008', lam_file_reses[index]]
#
#     for lam_file_res in iterlams:
#
#         for date in dates:
#
#             failed = 0
#
#             output_dir = '{}/reanalyzed/{}/'.format(root_dir, date)
#             if not os.path.exists(output_dir):
#                 os.makedirs(output_dir)
#
#             lamdate = date
#             critical_lambda_file = '{}/workdir/lambda/{}/lambdas_res={}.txt'.format(path_result, lamdate, lam_file_res)
#             critical_lambda_data = np.genfromtxt(critical_lambda_file)
#             xcrit = []
#             for i in range(122):
#                 xcrit = [i * 0.25 for i in range(122)]
#             ycrit = critical_lambda_data
#
#             if (resolution == '0.008' and date == '21_01_21_DNN1_024'):
#                 dataset = 100
#             else:
#                 dataset = datasets
#
#             for run_num in range(dataset):
#                 input_data_file = root_dir + date + '/sens_output_file_90CL_{}_resolution_{}.h5'.format(run_num, resolution)
#                 # Get the critical lambda data points, fit them to a spline curve.
#
#                 if not os.path.exists(input_data_file):
#                     failed += 1
#                     print('resolution = {}, iteration = {}, failed = {}'.format(resolution, run_num, failed))
#                     print(input_data_file)
#                     continue
#
#                 spline_xn = np.array([1., 3, 6, 10, 15, 21, 30])  # defines the locations of the knots
#                 # SplineFunc = LSQUnivariateSpline(critical_lambda_data[:,0],critical_lambda_data[:,1],t = spline_xn,k=3)
#                 SplineFunc = LSQUnivariateSpline(xcrit, ycrit, t=spline_xn, k=3)
#                 # Set some switches
#
#                 INCLUDE_EFFICIENCY_ERROR = False
#                 INCLUDE_BACKGROUND_SHAPE_ERROR = False
#                 PAR_LIMITS = True
#                 CONSTRAINTS = True
#                 DEBUG_PLOTTING = True
#
#                 input_df = pd.read_hdf(input_data_file)
#
#                 import time
#
#                 start_time = time.time()
#                 last_time = start_time
#
#                 new_limits = []
#
#                 for index, row in input_df.iterrows():
#                     # Next, find the hypothesis value for which lambda crosses the critical lambda curve
#                     if row['best_fit_covar']:
#                         xfit, yfit, crossing, crossing_idx = FindIntersectionByQuadraticInterpolation( \
#                             row['num_signal'][row['fixed_fit_acc_covar']], \
#                             row['lambda'][row['fixed_fit_acc_covar']], \
#                             SplineFunc)
#
#                         if DEBUG_PLOTTING:
#                             plt.plot(xfit, SplineFunc(xfit), '-b', label='90%CL spline')
#                             plt.plot(xcrit, ycrit, 'ob', linewidth=1, markersize=3, label='90%CL calculation')
#                             plt.plot(xfit, np.ones(len(xfit)) * 2.71, '--g', label='Wilks\'')
#                             plt.plot(xfit, yfit, '-r', label='Quadratic approx')
#                             plt.plot(row['num_signal'][row['fixed_fit_acc_covar']], \
#                                      row['lambda'][row['fixed_fit_acc_covar']], \
#                                      '-o', color='tab:blue', label='Actual fits')
#                             # plt.plot(xvals,lambdas,label='Actual fits')
#                             plt.plot(crossing, yfit[crossing_idx], 'ok')
#                             plt.xlim(0., 30.)
#                             plt.ylim(0., 5.)
#                             plt.xlabel('Num signal')
#                             plt.ylabel('Lambda')
#                             # plt.legend(loc='upper right')
#                             # plt.savefig('{}/example_fit_curve_DEBUGPLOT_{}.png'.format(output_dir, index), \
#                             #             dpi=200, bbox_inches='tight')
#                         new_limits.append(crossing)
#                     else:
#                         new_limits.append(-1.)
#
#                 print('\nDataset {} finished at {:4.4}s'.format(run_num, time.time() - last_time))
#                 print('Total time elapsed: {:4.4}s ({:4.4} min)\n\n\n'.format( \
#                     (time.time() - start_time), \
#                     (time.time() - start_time) / 60.))
#                 last_time = time.time()
#
#                 new_limits = np.array(new_limits)
#
#                 input_df['90CL_crossing_EXACT'] = new_limits
#
#                 ##########################################################################
#                 # Write the output dataframe to a file
#                 output_df = input_df
#
#                 input_data_file_name = input_data_file.split('/')[-1]
#
#                 # print(output_df.head())
#                 print('Saving file to output directory: {}'.format(output_dir))
#                 output_df.to_hdf('{}/reanalyzed_lres=_{}_{}'.format( \
#                     output_dir, lam_file_res, input_data_file_name), \
#                     key='df')
#
# print('Elapsed: {:4.4}s'.format(time.time() - start_time))


##########################################################################
##########################################################################

# input_data_file = sys.argv[1]
# critical_lambda_file = sys.argv[2]
# output_dir = sys.argv[3]

# Import useful libraries for analysis
import pandas as pd
import histlite as hl
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import LSQUnivariateSpline
from scipy.interpolate import UnivariateSpline

# Import the nEXO sensitivity classes
import nEXOFitWorkspace
import nEXOFitModel
import nEXOFitLikelihood

dnn = 'DNN1'
materialdb = "024"
dates = ["21_03_01_{}_{}".format(dnn, materialdb), "21_03_01_{}_023".format(dnn), "21_03_01_{}_025".format(dnn)]
materialdb = "024"
lamdate = "21_03_01_{}_{}".format(dnn, materialdb)
# resolutions = ['0.016', '0.017', '0.018']
# resolutions = ['0.011','0.015','0.018']
# resolutions = ['0.008']
resolutions = ['0.008', '0.009', '0.01', '0.011', '0.012', '0.013', '0.014', '0.015', '0.016', '0.017', '0.018']
datasets = 100

# path_home = '/Users/czyz1/lc-home/nexo_sensitivity/work'
path_result = '/p/lustre2/nexouser/czyz1'
# path_result = '/Users/czyz1/lc-nexouser'

root_dir = '{}/output/h5/'.format(path_result)
# lam_file_reses = ['0.011', '0.015', '0.018']
lam_file_reses = ['0.008']

# input_data_file = sys.argv[1]
# critical_lambda_file = sys.argv[2]
# output_dir = sys.argv[3]

for date in dates:

    if date.split('_')[-1] != '024':
        resolutions = ['0.008']

    iterlams = ['0.008']
    # else:
    #     iterlams = ['0.008', lam_file_reses[index]]

    for lam_file_res in iterlams:

        for index, resolution in enumerate(resolutions):

            failed = 0

            output_dir = '{}/reanalyzed/{}/'.format(root_dir, date)
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            critical_lambda_file = '{}/workdir/lambda/{}/lambdas_res={}.txt'.format(path_result, lamdate, lam_file_res)
            critical_lambda_data = np.genfromtxt(critical_lambda_file)
            xcrit = []
            for i in range(122):
                xcrit = [i * 0.25 for i in range(122)]
            ycrit = critical_lambda_data

            if (resolution == '0.008' and date == '21_01_21_DNN1_024'):
                dataset = 100
            else:
                dataset = datasets

            for run_num in range(dataset):
                input_data_file = root_dir + date + '/sens_output_file_90CL_{}_resolution_{}.h5'.format(run_num, resolution)
                # Get the critical lambda data points, fit them to a spline curve.

                if not os.path.exists(input_data_file):
                    failed += 1
                    print('resolution = {}, iteration = {}, failed = {}'.format(resolution, run_num, failed))
                    print(input_data_file)
                    continue

                # spline_xn = np.array([1., 3, 6, 10, 15, 21, 30])  # defines the locations of the knots
                spline_xn = np.array([.75, 2.25, 4.5, 10, 15, 21, 30])  # defines the locations of the knots
                # SplineFunc = LSQUnivariateSpline(critical_lambda_data[:,0],critical_lambda_data[:,1],t = spline_xn,k=3)
                SplineFunc = LSQUnivariateSpline(xcrit, ycrit, t=spline_xn, k=3)
                # Set some switches

                INCLUDE_EFFICIENCY_ERROR = False
                INCLUDE_BACKGROUND_SHAPE_ERROR = False
                PAR_LIMITS = True
                CONSTRAINTS = True
                DEBUG_PLOTTING = True

                input_df = pd.read_hdf(input_data_file)

                import time

                start_time = time.time()
                last_time = start_time

                new_limits = []

                for index, row in input_df.iterrows():
                    # Next, find the hypothesis value for which lambda crosses the critical lambda curve
                    if row['best_fit_covar']:
                        xfit, yfit, crossing, crossing_idx = FindIntersectionByQuadraticInterpolation( \
                            row['num_signal'][row['fixed_fit_acc_covar']], \
                            row['lambda'][row['fixed_fit_acc_covar']], \
                            SplineFunc)

                        if DEBUG_PLOTTING:
                            plt.plot(xfit, SplineFunc(xfit), '-b', label='90%CL spline')
                            plt.plot(xcrit, ycrit, 'ob', linewidth=1, markersize=3, label='90%CL calculation')
                            plt.plot(xfit, np.ones(len(xfit)) * 2.71, '--g', label='Wilks\'')
                            plt.plot(xfit, yfit, '-r', label='Quadratic approx')
                            plt.plot(row['num_signal'][row['fixed_fit_acc_covar']], \
                                     row['lambda'][row['fixed_fit_acc_covar']], \
                                     '-o', color='tab:blue', label='Actual fits')
                            # plt.plot(xvals,lambdas,label='Actual fits')
                            plt.plot(crossing, yfit[crossing_idx], 'ok')
                            plt.xlim(0., 30.)
                            plt.ylim(0., 5.)
                            plt.xlabel('Num signal')
                            plt.ylabel('Lambda')
                            # plt.legend(loc='upper right')
                            # plt.savefig('{}/example_fit_curve_DEBUGPLOT_{}.png'.format(output_dir, index), \
                            #             dpi=200, bbox_inches='tight')
                        new_limits.append(crossing)
                    else:
                        new_limits.append(-1.)

                print('\nDataset {} finished at {:4.4}s'.format(run_num, time.time() - last_time))
                print('Total time elapsed: {:4.4}s ({:4.4} min)\n\n\n'.format( \
                    (time.time() - start_time), \
                    (time.time() - start_time) / 60.))
                last_time = time.time()

                new_limits = np.array(new_limits)

                input_df['90CL_crossing_EXACT'] = new_limits

                ##########################################################################
                # Write the output dataframe to a file
                output_df = input_df

                input_data_file_name = input_data_file.split('/')[-1]

                # print(output_df.head())
                print('Saving file to output directory: {}'.format(output_dir))
                output_df.to_hdf('{}/{}'.format(output_dir, input_data_file_name), key='df')

print('Elapsed: {:4.4}s'.format(time.time() - start_time))



