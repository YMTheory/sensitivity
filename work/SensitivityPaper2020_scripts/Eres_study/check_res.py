# Import useful libraries for analysis
import sys
import os
import pandas as pd
import histlite as hl
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit

sys.path.append('../../modules')

# Import the nEXO sensitivity classes
import nEXOFitWorkspace
import nEXOFitModel
import nEXOFitLikelihood


def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def main():
    path_home = '/Users/czyz1/lc-home/nexo_sensitivity/work'
    path_result = '/Users/czyz1/lc-nexouser'

    config_loc = "{}/config/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1.yaml".format(path_home)
    date = '21_03_01'
    database_num = '023'
    # Set some switches
    lt_years = 10
    fig = None

    for resolution in ['0.008', '0.009', '0.01', '0.011', '0.012', '0.013', '0.014', '0.015', '0.016', '0.017', '0.018']:

        comp_loc = "{}/workdir/components_tables/{}/ComponentsTable_D-{}_Energy_Res={}.h5".\
            format(path_result, date, database_num, resolution)

        # Create the workspace
        workspace = nEXOFitWorkspace.nEXOFitWorkspace(config_loc)
        workspace.SetHandlingOfRadioassayData(fluctuate=True)

        workspace.LoadComponentsTableFromFile(comp_loc)
        workspace.livetime = lt_years * 365.25 * 24. * 60. * 60.
        workspace.CreateGroupedPDFs()

        # Create the likelihood object
        likelihood = nEXOFitLikelihood.nEXOFitLikelihood()
        likelihood.AddPDFDataframeToModel(workspace.df_group_pdfs, workspace.histogram_axis_names)

        initial_guess = likelihood.GetVariableValues()
        print(initial_guess)
        likelihood.model.UpdateVariables(initial_guess)
        likelihood.model.GenerateModelDistribution()
        likelihood.AddDataset(likelihood.model.GenerateDataset())

        ss_cut_dict = {'DNN': (0.15, 1), 'Energy (keV)': (700., 3500.), 'Standoff (mm)': (0., 650.)}
        ms_cut_dict = {'DNN': (0., 0.15), 'Energy (keV)': (700., 3500.), 'Standoff (mm)': (0., 650.)}

        # Set up the plotting parameters
        plt.rcParams.update({'font.size': 14})

        if not fig:
            fig, ax = plt.subplots(1, 1, figsize=(12, 10))

        weight = 1000000  # var['Value']
        print('name = ' + 'Num_FullLXeBb0n' + ', weight = ' + str(weight))
        # ['Rn222_FieldRingRadon', 'Rn222_CathodeRadon', 'Rn222_ActiveLXe', 'Rn222_InactiveLXe']
        ss_pdf = likelihood.model.GetSlicedDistribution(ss_cut_dict, var_name='Num_FullLXeBb0n', verbose=False)
        ms_pdf = likelihood.model.GetSlicedDistribution(ms_cut_dict, var_name='Num_FullLXeBb0n', verbose=False)

        print('Plotting {}'.format('Num_FullLXeBb0n'))

        ss_sum = hl.hist([np.array([0.]), np.array([0.]), np.array([0.])], \
                         bins=ss_pdf.bins)
        ms_sum = hl.hist([np.array([0.]), np.array([0.]), np.array([0.])], \
                         bins=ms_pdf.bins)

        hn = (weight * ss_pdf).project([1])
        # hl.plot1d(ax, (weight * ss_pdf).project([1]), label='Resolution = {}'.format(resolution))
        # hl.plot1d(likelihood.ax[0, 1], (weight * ms_pdf).project([1]))
        # hl.plot1d(likelihood.ax[1, 1], (weight * ms_pdf).project([2]))
        # hl.plot1d(likelihood.ax[1, 0], (weight * ss_pdf).project([2]))


        coeff, var_matrix = curve_fit(gauss, hn.centers[0], hn.values, p0=[1, 2458, 20])
        hist_fit = gauss(hn.centers[0], *coeff)
        # plt.plot(hn.centers[0], hist_fit, label='Res_fit = {}'.format(coeff[2]/2458))
        plt.plot(hn.centers[0], hist_fit, label='Resolution = {}'.format(resolution))

        # X = np.linspace(2300, 2600, 2600-2300)
        #
        # params, cov = hn.curve_fit(lambda x=X, loc=2458, scale=1: stats.norm.pdf(x, loc, scale))
        # ax.plot(X, stats.norm.pdf(X, *params))

        fig.legend(ncol=4, facecolor=(1., 1., 1.), framealpha=1., loc='upper left')
        # likelihood.ax[0,0].set_ylim(1e-2,1e7)
        ax.set_xlim(2300., 2620.)
        ax.set_ylabel('Counts')
        ax.set_xlabel('Energy (keV)')
        # likelihood.ax[0,0].set_yscale('log')
        # likelihood.ax[0, 1].set_ylim(1e-2, 1e7)
        # likelihood.ax[0, 1].set_xlim(700., 3500.)
        # likelihood.ax[0, 1].set_ylabel('Counts')
        # likelihood.ax[0, 1].set_xlabel('Energy (keV)')
        # likelihood.ax[0, 1].set_yscale('log')
        # likelihood.ax[1, 1].set_ylim(1e-2, 1e7)
        # likelihood.ax[1, 1].set_xlim(0., 640.)
        # likelihood.ax[1, 1].set_ylabel('Counts')
        # likelihood.ax[1, 1].set_xlabel('Standoff (mm)')
        # likelihood.ax[1, 1].set_yscale('log')
        # likelihood.ax[1, 0].set_ylim(1e-2, 1e7)
        # likelihood.ax[1, 0].set_xlim(0., 640.)
        # likelihood.ax[1, 0].set_ylabel('Counts')
        # likelihood.ax[1, 0].set_xlabel('Standoff (mm)')
        # likelihood.ax[1, 0].set_yscale('log')


    plt.show()
    plt.savefig('{}/output/plots/{}/pdfs_{}.png'.format(path_result, date, resolution), dpi=200, bbox_inches='tight')


if __name__ == "__main__":
    main()