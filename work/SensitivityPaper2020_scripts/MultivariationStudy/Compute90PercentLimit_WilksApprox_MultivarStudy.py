# Import sys, then tell python where to find the nEXO-specific classes
import sys
sys.path.append('../../../modules')

# Import useful libraries
import argparse
import pathlib
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import time
import hashlib

# Import the nEXO sensitivity classes
import nEXOFitWorkspace
import nEXOFitLikelihood


def get_parser():
    parser = argparse.ArgumentParser(description='Compute 90% limit with Wilks Approximation')
    parser.add_argument("job_id_num", type=int, help='job ID number')
    parser.add_argument("--bkg_shape_err", type=float, default=0., help='background shape error parameter')
    parser.add_argument("-n", "--num_datasets", type=int, default=1, help='Number of datasets to generate')
    parser.add_argument("input_table", type=pathlib.Path, help='Input table')
    parser.add_argument("output_dir", type=pathlib.Path, help='Output directory')
    parser.add_argument("-e", "--energy_res", type=float, default=0.008, help='Energy resolution')
    parser.add_argument("-d", "--dnn_scale_factor", type=float, default=0., help='DNN distribution smearing factor')
    parser.add_argument("-b", "--bkg_scale_factor", type=float, default=1.,
                        help='Gamma background scale factor (Rn222 excluded)')
    parser.add_argument("-x", "--xe137_scale_factor", type=float, default=1., help='Xe-137/Ar-42 scale factor')
    parser.add_argument("-r", "--rn222_scale_factor", type=float, default=1., help='Rn-222 scale factor')
    parser.add_argument("-c", "--workspace_config", type=pathlib.Path,
                        default='/p/vast1/nexo/sensitivity2020/pdfs/config_files/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1.yaml',
                        help="nEXOFitWorkspace configuration file")
    parser.add_argument("--debug", dest='debug', action='store_true', help='Turn on debug plotting')
    parser.set_defaults(debug=False)
    return parser


def FindIntersectionByQuadraticInterpolationWilks(xvals, yvals):
    if len(xvals) != len(yvals):
        print('ERROR: need same length arrays for ' +
              'quadratic interpolation')
        raise ValueError
    print('Xvals:')
    print(xvals)
    print('Yvals:')
    print(yvals)

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
            print(f'Quadratic fit: {p}')
            yfit = p[0] * xfit ** 2 + p[1] * xfit + p[2]
        except np.RankWarning:
            p = np.polyfit(xvals[mask], yvals[mask], 1)
            print(f'Linear fit: {p}')
            yfit = p[0] * xfit + p[1]
        ythreshold = np.ones(len(yfit)) * 2.706  # This is the Wilks' approx.
        crossing_idx = np.where((yfit - ythreshold) > 0.)[0][0]
        crossing = xfit[crossing_idx]
    else:
        yfit = np.zeros(len(xfit))
        crossing_idx = -1
        crossing = -1.

    return xfit, yfit, crossing, crossing_idx


if __name__ == "__main__":
    arg_parser = get_parser()
    args = arg_parser.parse_args()
    for arg in vars(args):
        print(arg, getattr(args, arg))

    # Set some switches
    INCLUDE_EFFICIENCY_ERROR = False
    INCLUDE_BACKGROUND_SHAPE_ERROR = False
    PAR_LIMITS = True
    CONSTRAINTS = True

    # Create the workspace and load smeared energy resolution components table
    workspace = nEXOFitWorkspace.nEXOFitWorkspace(args.workspace_config)
    workspace.LoadComponentsTableFromFile(args.input_table)
    # workspace.LoadComponentsTableFromFile(args.input_table + args.energy_res + '.h5')

    # from tabulate import tabulate
    # print(tabulate(workspace.df_components, headers='keys', tablefmt='psql'))

    for index, row in workspace.df_components.iterrows():
        # Scale the gamma ray background components, except radon
        isotopes_to_leave_alone = ['Ar42', 'Xe137', 'bb2n', 'bb0n', 'B8nu', 'Rn222', ]  # just for bookkeeping
        isotopes_to_scale = ['K40', 'Co60', 'Al26', 'Th232', 'U238', 'Cs137']
        # The format is <isotope>_<part>, e.g. "Th232_HVCables"
        if row['PDFName'].split('_')[0] in isotopes_to_scale:
            print(f'Scaling {row["PDFName"]}...')
            workspace.df_components.loc[index, 'SpecActiv'] = args.bkg_scale_factor * row['SpecActiv']
            workspace.df_components.loc[index, 'SpecActivErr'] = args.bkg_scale_factor * row['SpecActivErr']

        # Scale the Xe137 and Ar42 components.
        if 'Xe137' in row['PDFName'] or 'Ar42' in row['PDFName']:
            print(f'Scaling {row["PDFName"]}...')
            workspace.df_components.loc[index, 'SpecActiv'] = args.xe137_scale_factor * row['SpecActiv']
            workspace.df_components.loc[index, 'SpecActivErr'] = args.xe137_scale_factor * row['SpecActivErr']

    workspace.CreateGroupedPDFs()

    # Define the ROI within the workspace
    roi_dict = {'DNN': [0.85, 1.],
                'Energy (keV)': [2434., 2480.],
                'Standoff (mm)': [104.5, 650.]}
    workspace.DefineROI(roi_dict)

    # Create the likelihood object
    likelihood = nEXOFitLikelihood.nEXOFitLikelihood()
    likelihood.AddPDFDataframeToModel(workspace.df_group_pdfs, workspace.histogram_axis_names)

    if INCLUDE_EFFICIENCY_ERROR:
        likelihood.model.IncludeSignalEfficiencyVariableInFit(True)
    if INCLUDE_BACKGROUND_SHAPE_ERROR:
        likelihood.model.IncludeBackgroundShapeVariableInFit(True)

    initial_guess = likelihood.GetVariableValues()

    # Scale the Rn222 component according to the input value
    rn222_idx = likelihood.model.GetVariableIndexByName('Rn222')
    initial_guess[rn222_idx] *= args.rn222_scale_factor

    # Update the model in the likelihood object
    likelihood.model.UpdateVariables(initial_guess)
    likelihood.model.GenerateModelDistribution()

    # Print out the number of events in the ROI
    total_bkg_in_roi = likelihood.model.GetIntegralInBinRange(workspace.GetROIBinIndices())
    print('\n****************************************************************************************')
    print('Variable list after changes to components and resolution:')
    likelihood.PrintVariableList()
    print(f'\n\nTotal ROI background: '
          f'{likelihood.model.GetIntegralInBinRange(workspace.GetROIBinIndices()):4.4} events\n')
    print('Contribution from each component in ROI:')
    for component in likelihood.model.variable_list:
        if 'Shape' in component['Name']:
            continue
        num_counts_in_roi = likelihood.model.GetComponentIntegralInBinRange(
            component['Name'], workspace.GetROIBinIndices())
        print(f'{component["Name"] + ":":<20}\t'
              f'{num_counts_in_roi:>10.4}\t'
              f'{int(1000 * num_counts_in_roi / total_bkg_in_roi) / 10.:>10.4}%')
    print('****************************************************************************************\n')

    # Add an initial dataset
    likelihood.AddDataset(likelihood.model.GenerateDataset())

    # Set limits so that none of the PDFs can go negative in the fit.
    # (except the signal PDF)
    if PAR_LIMITS:
        for var in likelihood.model.variable_list:
            if 'Bb0n' in var['Name']:
                likelihood.SetVariableLimits(var['Name'],
                                             lower_limit=-15.,
                                             upper_limit=100.)
            elif 'Co60' in var['Name']:
                likelihood.SetVariableLimits(var['Name'],
                                             lower_limit=0.,
                                             upper_limit=var['Value'] * 10.)
            else:
                likelihood.SetVariableLimits(var['Name'],
                                             lower_limit=var['Value'] * 0.1,
                                             upper_limit=var['Value'] * 10.)

    likelihood.SetFractionalMinuitInputError('Num_FullLXeBb0n', 0.01 / 0.0001)

    ##########################################################################
    # Here's where the calculation loop begins.
    workspace.SetHandlingOfRadioassayData(fluctuate=True)

    num_hypotheses = 25
    xvals = np.array([])

    lambdas = np.zeros(num_hypotheses)
    num_iterations = np.zeros(num_hypotheses)
    best_fit_converged = True
    crossing = -1
    output_df_list = []

    start_time = time.time()
    last_time = start_time

    for j in range(0, args.num_datasets):

        # Initialize (or reset) all my output variables.
        converged = True
        num_iterations = np.ones(num_hypotheses)
        lambdas = np.zeros(num_hypotheses)
        xvals = np.zeros(num_hypotheses)
        fixed_fit_converged = np.array([], dtype=bool)
        fixed_fit_covar = np.array([], dtype=bool)
        crossing = -1
        output_row = dict()

        # Redo the grouping, which fluctuates the radioassay values within their uncertainties.
        workspace.CreateGroupedPDFs()
        likelihood.AddPDFDataframeToModel(workspace.df_group_pdfs,
                                          workspace.histogram_axis_names,
                                          replace_existing_variables=False)

        rn222_idx = likelihood.model.GetVariableIndexByName('Rn222')
        likelihood.model.variable_list[rn222_idx]['Value'] *= args.rn222_scale_factor

        initial_guess = likelihood.GetVariableValues()
        likelihood.model.GenerateModelDistribution()
        likelihood.AddDataset(likelihood.model.GenerateDataset())

        likelihood.SetAllVariablesFloating()

        # Fix the Co60 parameter
        # likelihood.SetVariableFixStatus('Num_FullTPC_Co60',True)

        if CONSTRAINTS:
            rn222_idx = likelihood.GetVariableIndex('Rn222')
            # Fluctuate Rn222 constraint
            rn222_constraint_val = (np.random.randn() * 0.1 + 1) * initial_guess[rn222_idx]
            # Set Rn222 constraint
            likelihood.SetGaussianConstraintAbsolute(likelihood.model.variable_list[rn222_idx]['Name'],
                                                     rn222_constraint_val,
                                                     0.1 * initial_guess[rn222_idx])
            b8_idx = likelihood.GetVariableIndex('B8')
            # Fluctuate B8nu constraint
            b8_constraint_val = (np.random.randn() * 0.1 + 1) * initial_guess[b8_idx]
            # Set B8nu constraint
            likelihood.SetGaussianConstraintAbsolute(likelihood.model.variable_list[b8_idx]['Name'],
                                                     b8_constraint_val,
                                                     0.1 * initial_guess[b8_idx])

        if INCLUDE_EFFICIENCY_ERROR:
            eff_idx = likelihood.GetVariableIndex('Signal_Efficiency')
            eff_constraint_val = (np.random.randn() * eff_err + 1) * initial_guess[eff_idx]
            likelihood.SetGaussianConstraintAbsolute(likelihood.model.variable_list[eff_idx]['Name'],
                                                     eff_constraint_val,
                                                     eff_err)
        if INCLUDE_BACKGROUND_SHAPE_ERROR:
            bkg_shape_idx = likelihood.GetVariableIndex('Background_Shape_Error')
            bkg_shape_constraint_val = np.random.randn() * args.bkg_shape_err
            likelihood.SetGaussianConstraintAbsolute(likelihood.model.variable_list[bkg_shape_idx]['Name'],
                                                     bkg_shape_constraint_val,
                                                     args.bkg_shape_err)

        print(f'\n\nRunning dataset {j}...\n')
        likelihood.PrintVariableList()

        print('\nConstraints:')
        for constraint in likelihood.model.constraints:
            print(f'\t{constraint}')
        print('\n')

        # Fix the signal parameter to different values, and calculate lambda
        for i in (range(0, num_hypotheses)):

            # Fix the 0nu parameter to a specific hypothesis
            if args.xe137_scale_factor > 9. or args.bkg_shape_err > 9.:
                signal_hypothesis = float(i) * 2.2 + 0.000001
            else:
                signal_hypothesis = float(i) * 1.4 + 0.000001
            signal_idx = likelihood.GetVariableIndex('Bb0n')
            initial_values = np.copy(initial_guess)
            initial_values[signal_idx] = signal_hypothesis
            xvals[i] = signal_hypothesis

            ###########################################################################
            # All the exciting stuff happens here!
            lambda_fit_result = likelihood.ComputeLambdaForPositiveSignal(
                initial_values=initial_values,
                signal_name='Bb0n',
                signal_expectation=0.,
                print_level=1,
                fixed_fit_signal_value=signal_hypothesis)
            ###########################################################################

            # Store all the important quantities
            lambdas[i] = lambda_fit_result['lambda']
            fixed_fit_converged = np.append(fixed_fit_converged, lambda_fit_result['fixed_fit_converged'])
            fixed_fit_covar = np.append(fixed_fit_covar, lambda_fit_result['fixed_fit_covar'])
            num_iterations[i] = lambda_fit_result['fixed_fit_iterations']

            print('After fit ends:')
            likelihood.PrintVariableList()
            print('\n')

        # Next, find the hypothesis value for which lambda crosses the critical lambda curve
        if lambda_fit_result['best_fit_covar']:
            xfit, yfit, crossing, crossing_idx = FindIntersectionByQuadraticInterpolationWilks(
                xvals[fixed_fit_covar],
                lambdas[fixed_fit_covar])

            if args.debug:
                plt.clf()
                plt.plot(xfit, np.ones(len(xfit)) * 2.706, '-b', label='90%CL spline')
                plt.plot(xfit, yfit, '-r', label='Quadratic approx')
                plt.plot(xvals[fixed_fit_converged],
                         lambdas[fixed_fit_converged],
                         '-o', color='tab:blue', label='Actual fits')
                # plt.plot(xvals,lambdas,label='Actual fits')
                plt.plot(crossing, yfit[crossing_idx], 'ok')
                plt.xlim(0., 30.)
                plt.ylim(0., 7.)
                plt.xlabel('Num signal')
                plt.ylabel('Lambda')
                plt.legend(loc='upper right')
                plt.savefig(f'{args.output_dir}/example_fit_curve_DEBUGPLOT_{j}.png',
                            dpi=200, bbox_inches='tight')

        output_row['num_signal'] = xvals
        output_row['lambda'] = lambdas
        output_row['fixed_fit_converged'] = fixed_fit_converged
        output_row['fixed_fit_acc_covar'] = fixed_fit_covar
        output_row['90CL_crossing'] = crossing
        output_row['num_iterations'] = num_iterations
        # The "best_fit" quantities in lambda_fit_result should be the same for
        # every lambda calculation, so it's okay if we use the most recent one
        output_row['best_fit_converged'] = lambda_fit_result['best_fit_converged']
        output_row['best_fit_covar'] = lambda_fit_result['best_fit_covar']
        output_row['best_fit_iterations'] = lambda_fit_result['best_fit_iterations']
        output_row['best_fit_parameters'] = lambda_fit_result['best_fit_parameters']
        output_row['best_fit_errors'] = lambda_fit_result['best_fit_errors']
        output_row['best_fit_nll'] = lambda_fit_result['best_fit_nll']
        output_row['fixed_fit_parameters'] = lambda_fit_result['fixed_fit_parameters']
        output_row['fixed_fit_errors'] = lambda_fit_result['fixed_fit_errors']
        output_row['input_parameters'] = initial_guess
        # output_row['dataset'] = likelihood.dataset

        output_df_list.append(output_row)

        print(f'\nDataset {j} finished at {time.time() - last_time:4.4}s')
        print(f'Total time elapsed: {(time.time() - start_time):4.4}s '
              f'({(time.time() - start_time) / 60.:4.4} min)\n\n\n')
        last_time = time.time()

    # Write the output dataframe and metadata to files
    output_df = pd.DataFrame(output_df_list)
    # print(output_df.head())

    s = f'Xe137:{args.xe137_scale_factor:0>4.4f} ' + \
        f'Rn222:{args.rn222_scale_factor:0>4.4f} ' + \
        f'DNN:{args.dnn_scale_factor:0>4.4f} ' + \
        f'Bkg:{args.bkg_scale_factor:0>4.4f} ' + \
        f'ERes:{args.energy_res:0>4.4f}'
    s_hash = hashlib.md5(s.encode('utf-8')).hexdigest()[:6].upper()

    print(f'Saving file to output directory: {args.output_dir}')
    output_df.to_hdf(f'{args.output_dir}/'
                     f'sens_output_file_multivarstudy_{s_hash}'
                     f'x_90CL_{args.job_id_num:03}.h5',
                     key='df')

    with open(f'{args.output_dir}/{s_hash}.txt', 'a') as file:
        print(f'{args.job_id_num:03}', s, file=file)
        print(f'{args.job_id_num:03}', ' '.join(f'{k}={v}' for k, v in vars(args).items()), file=file)

    print(f'Elapsed: {time.time() - start_time:4.4}s')
