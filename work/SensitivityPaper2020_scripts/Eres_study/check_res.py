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

config_loc = "/p/lustre2/czyz1/nexo_sensitivity/work/config/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1_v9wAr42.yaml"

# Set some switches
lt_years = 10

for resolution in ['0.008', '0.009', '0.01', '0.011', '0.012', '0.013', '0.014', '0.015', '0.016']:

    comp_loc = "/p/lustre2/nexouser/czyz1/workdir/components_tables/ComponentsTable_D-024_wAr42_Energy_Res={}.h5".\
        format(resolution)

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

    ss_cut_dict = {'DNN': (0.2, 1), 'Energy (keV)': (700., 3500.), 'Standoff (mm)': (0., 650.)}
    ms_cut_dict = {'DNN': (0., 0.2), 'Energy (keV)': (700., 3500.), 'Standoff (mm)': (0., 650.)}

    # Set up the plotting parameters
    plt.rcParams.update({'font.size': 14})

    likelihood.fig, likelihood.ax = plt.subplots(2, 2, figsize=(12, 10))

    weight = 1000000  # var['Value']
    print('name = ' + 'Num_FullLXeBb0n' + ', weight = ' + str(weight))
    ss_pdf = likelihood.model.GetSlicedDistribution(ss_cut_dict, var_name='Num_FullLXeBb0n', verbose=False)
    ms_pdf = likelihood.model.GetSlicedDistribution(ms_cut_dict, var_name='Num_FullLXeBb0n', verbose=False)

    print('Plotting {}'.format('Num_FullLXeBb0n'))

    ss_sum = hl.hist([np.array([0.]), np.array([0.]), np.array([0.])], \
                     bins=ss_pdf.bins)
    ms_sum = hl.hist([np.array([0.]), np.array([0.]), np.array([0.])], \
                     bins=ms_pdf.bins)

    hn = (weight * ss_pdf).project([1])
    hl.plot1d(likelihood.ax[0, 0], (weight * ss_pdf).project([1]), label='Num_FullLXeBb0n')
    hl.plot1d(likelihood.ax[0, 1], (weight * ms_pdf).project([1]))
    hl.plot1d(likelihood.ax[1, 1], (weight * ms_pdf).project([2]))
    hl.plot1d(likelihood.ax[1, 0], (weight * ss_pdf).project([2]))
    X = np.linspace(2300, 2600, 2600-2300)

    params, cov = hn.curve_fit(lambda x=X, loc=2458, scale=1: stats.norm.pdf(x, loc, scale))
    likelihood.ax[0, 0].plot(X, stats.norm.pdf(X, *params))

    likelihood.fig.legend(ncol=4, facecolor=(1., 1., 1.), framealpha=1., loc='upper center')
    # likelihood.ax[0,0].set_ylim(1e-2,1e7)
    likelihood.ax[0, 0].set_xlim(2300., 2600.)
    likelihood.ax[0, 0].set_ylabel('Counts')
    likelihood.ax[0, 0].set_xlabel('Energy (keV)')
    # likelihood.ax[0,0].set_yscale('log')
    likelihood.ax[0, 1].set_ylim(1e-2, 1e7)
    likelihood.ax[0, 1].set_xlim(700., 3500.)
    likelihood.ax[0, 1].set_ylabel('Counts')
    likelihood.ax[0, 1].set_xlabel('Energy (keV)')
    likelihood.ax[0, 1].set_yscale('log')
    likelihood.ax[1, 1].set_ylim(1e-2, 1e7)
    likelihood.ax[1, 1].set_xlim(0., 640.)
    likelihood.ax[1, 1].set_ylabel('Counts')
    likelihood.ax[1, 1].set_xlabel('Standoff (mm)')
    likelihood.ax[1, 1].set_yscale('log')
    likelihood.ax[1, 0].set_ylim(1e-2, 1e7)
    likelihood.ax[1, 0].set_xlim(0., 640.)
    likelihood.ax[1, 0].set_ylabel('Counts')
    likelihood.ax[1, 0].set_xlabel('Standoff (mm)')
    likelihood.ax[1, 0].set_yscale('log')

    plt.savefig('/p/lustre2/nexouser/czyz1/output/plots/pdfs_{}.png'.format(resolution), dpi=200, bbox_inches='tight')

    plt.show()
