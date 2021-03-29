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
dates = ["21_03_11_{}_{}".format(dnn, materialdb)]
lamdate = "21_03_22_{}_{}".format(dnn, materialdb)
livetimes = ['0.5', '1.0', '2.0', '5.0', '10.0']
datasets = 100

# path_home = '/Users/czyz1/lc-home/nexo_sensitivity/work'
path_result = '/p/lustre2/nexouser/czyz1'
# path_result = '/Users/czyz1/lc-nexouser'

root_dir = '{}/output/h5/'.format(path_result)
lam_file_res = '0.008'

# input_data_file = sys.argv[1]
# critical_lambda_file = sys.argv[2]
# output_dir = sys.argv[3]

for index, lt in enumerate(livetimes):

    resolution = '0.008'

    for date in dates:

        failed = 0

        output_dir = '{}/reanalyzed/{}/'.format(root_dir, date)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        critical_lambda_file = '{}/workdir/lambda/{}/lambdas_lt={}.txt'.format(path_result, lamdate, lt)
        critical_lambda_data = np.genfromtxt(critical_lambda_file)
        xcrit = []
        for i in range(122):
            xcrit = [i * 0.25 for i in range(122)]
        ycrit = critical_lambda_data

        for run_num in range(datasets):
            input_data_file = root_dir + date + '/sens_output_file_90CL_{}_livetime_{}.h5'.format(run_num, lt)
            # Get the critical lambda data points, fit them to a spline curve.

            if not os.path.exists(input_data_file):
                failed += 1
                print('resolution = {}, iteration = {}, failed = {}'.format(resolution, run_num, failed))
                print(input_data_file)
                continue

            spline_xn = np.array([1., 3, 6, 10, 15, 21, 30])  # defines the locations of the knots
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
            output_df.to_hdf('{}/{}'.format( \
                output_dir, input_data_file_name), key='df')

print('Elapsed: {:4.4}s'.format(time.time() - start_time))



