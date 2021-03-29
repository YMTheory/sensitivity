# Import sys, then tell python where to find the nEXO-specific classes
import sys
import os

sys.path.append('../../modules')

######################################################################
# Load the arguments:
if len(sys.argv) != 9:
    print('\nERROR: incorrect number of arguments.\n')
    print('Usage:')
    print('\tpython Compute90PercentLimit_Wilks_Cluster.py ' + \
        '<iteration_num> <bkg_shape_err_parameter> <num_datasets_to_generate> <output_directory>\n\n')
    sys.exit()

iteration_num = int(sys.argv[1])
bkg_shape_err = float(sys.argv[2])
num_datasets = int(sys.argv[3])
output_dir = sys.argv[4]
num_hypotheses = int(sys.argv[5])
date = sys.argv[6]
config_loc = sys.argv[7]
comp_loc = sys.argv[8]
##########################################################################
def FindIntersectionByQuadraticInterpolationWilks( xvals, yvals ):

	if len(xvals)!=len(yvals):
		print('ERROR: need same length arrays for ' +\
			'quadratic interpolation')
		raise ValueError
	print('Xvals:')
	print(xvals)
	print('Yvals:')
	print(yvals)

	# First, select only lambda values on the upward slope, so we find
	# the upper limit
	mask = np.zeros(len(yvals),dtype=bool)
	mask[1:] = ( yvals[1:] - yvals[:-1] ) > 0.
	# Next, select only values near the critical lambda threshold (~2.7)
	mask = mask&(yvals>0.5)&(yvals<6.)

	xfit = np.linspace(0.,40.,800)

	if len( xvals[mask] ) > 0:
		try:
			p = np.polyfit( xvals[mask], yvals[mask], 2)
			print('Quadratic fit: {}'.format(p))
			yfit = p[0]*xfit**2 + p[1]*xfit + p[2]
		except np.RankWarning:
			p = np.polyfit( xvals[mask], yvals[mask], 1)
			print('Linear fit: {}'.format(p))
			yfit = p[0]*xfit + p[1]
		ythreshold = np.ones(len(yfit))*2.706 # This is the Wilks' approx.
		crossing_idx = np.where( (yfit - ythreshold)>0. )[0][0]
		crossing = xfit[crossing_idx]
	else:
		yfit = np.zeros(len(xfit))
		crossing_idx = -1
		crossing = -1.

	return xfit, yfit, crossing, crossing_idx
##########################################################################
##########################################################################


# Import useful libraries for analysis
import pandas as pd
import histlite as hl
import numpy as np
from matplotlib import pyplot as plt


# Import the nEXO sensitivity classes
import nEXOFitWorkspace
import nEXOFitModel
import nEXOFitLikelihood



# Set some switches
INCLUDE_EFFICIENCY_ERROR = False
INCLUDE_BACKGROUND_SHAPE_ERROR = False
PAR_LIMITS = True
DEBUG_PLOTTING = False
lt_years = 10

for resolution in ['0.008']:
	# Create the workspace
	workspace = nEXOFitWorkspace.nEXOFitWorkspace(config_loc)
	workspace.SetHandlingOfRadioassayData(fluctuate=True)
	
	workspace.ConnectToMaterialsDB()
	workspace.materialsDB.PrintAllGeometryTags()
	# workspace.CreateComponentsTableFromMaterialsDB('D-024',
        #                                                histograms_file='/p/lustre2/czyz1/nexo_sensitivity/work/histogram_files/Baseline2019_Histograms_D-024_Energy_Res={}.h5'.format(resolution))
	#                                               histograms_file='/p/vast1/nexo/sensitivity2020/pdfs/component_tables/ComponentsTable_D-024_Optimized_DNN_Standoff_Binning_version1.h5')
	# print(workspace.df_components.head())
	workspace.LoadComponentsTableFromFile(comp_loc+resolution+'.h5')
	workspace.livetime = lt_years * 365.25 * 24. * 60. * 60.
	workspace.CreateGroupedPDFs()


	# Create the likelihood object
	likelihood = nEXOFitLikelihood.nEXOFitLikelihood()
	likelihood.AddPDFDataframeToModel(workspace.df_group_pdfs, workspace.histogram_axis_names )


	if INCLUDE_EFFICIENCY_ERROR:
		likelihood.model.IncludeSignalEfficiencyVariableInFit(True)
	if INCLUDE_BACKGROUND_SHAPE_ERROR:
		likelihood.model.IncludeBackgroundShapeVariableInFit(True)


	initial_guess = likelihood.GetVariableValues()
	print(initial_guess)
	likelihood.model.UpdateVariables(initial_guess)
	likelihood.model.GenerateModelDistribution()
	likelihood.AddDataset( likelihood.model.GenerateDataset() )


	# Set limits so that none of the PDFs go negative in the fit.
	# (except the signal PDF)
	if PAR_LIMITS:
		for var in likelihood.model.variable_list:
			if 'Bb0n' in var['Name']:
				likelihood.SetVariableLimits( var['Name'], \
										  lower_limit = -15., \
										  upper_limit = 100.)
			elif 'Background_Shape_Error' in var['Name']:
				likelihood.SetVariableLimits( var['Name'], \
										  lower_limit = -100., \
										  upper_limit = 100.)
			else:
				likelihood.SetVariableLimits( var['Name'], \
										  lower_limit = 0., \
										  upper_limit = var['Value']*10.)

	likelihood.SetFractionalMinuitInputError('Num_FullLXeBb0n', 0.01/0.0001)


	##########################################################################
	# Here's where the calculation loop begins.

	# num_hypotheses = 15
	xvals = np.array([])

	lambdas = np.zeros(num_hypotheses)
	num_iterations = np.zeros(num_hypotheses)
	best_fit_converged = True
	crossing = -1
	output_df_list = []

	import time

	start_time = time.time()
	last_time = start_time

	if DEBUG_PLOTTING:
		ss_cut_dict = {'DNN': (0.85, 1), 'Energy (keV)': (700., 3500.), 'Standoff (mm)': (0., 650.)}
		ms_cut_dict = {'DNN': (0., 0.85), 'Energy (keV)': (700., 3500.), 'Standoff (mm)': (0., 650.)}

		likelihood.PlotModelDistributions(ss_cut_dict, ms_cut_dict, \
									  output_filename='/p/lustre2/nexouser/czyz1/output/plots/'+date+'_test.png', \
									  save=True, show=True, plot_data=True)

	for j in range(0,num_datasets):

		# Initialize (or reset) all my output variables.
		converged = True
		num_iterations = np.ones(num_hypotheses)
		lambdas = np.zeros(num_hypotheses)
		xvals = np.zeros(num_hypotheses)
		fixed_fit_converged = np.array([],dtype=bool)
		fixed_fit_covar = np.array([],dtype=bool)
		crossing = -1
		output_row = dict()

		workspace.CreateGroupedPDFs()
		likelihood.AddPDFDataframeToModel(workspace.df_group_pdfs, workspace.histogram_axis_names,
										  replace_existing_variables=False)

		# likelihood.model.UpdateVariables(initial_guess)
		likelihood.model.GenerateModelDistribution()
		likelihood.AddDataset( likelihood.model.GenerateDataset() )

		likelihood.SetAllVariablesFloating()

		# Fix the Co60 parameter
		# likelihood.SetVariableFixStatus('Num_FullTPC_Co60',True)

		RN_CONSTRAINTS=True
		if RN_CONSTRAINTS:
			rn222_idx = likelihood.GetVariableIndex('Rn222')
			# Fluctuate Rn222 constraint
			rn222_constraint_val = (np.random.randn()*0.1 + 1)*initial_guess[rn222_idx]
			# Set Rn222 constraint
			likelihood.SetGaussianConstraintAbsolute(likelihood.model.variable_list[rn222_idx]['Name'],\
								 rn222_constraint_val, \
													 0.1 * initial_guess[rn222_idx])

		if INCLUDE_EFFICIENCY_ERROR:
					eff_idx = likelihood.GetVariableIndex('Signal_Efficiency')
					eff_constraint_val = (np.random.randn()*eff_err + 1)* initial_guess[eff_idx]
					likelihood.SetGaussianConstraintAbsolute(likelihood.model.variable_list[eff_idx]['Name'],\
															eff_constraint_val,\
															eff_err)
		if INCLUDE_BACKGROUND_SHAPE_ERROR:
			bkg_shape_idx = likelihood.GetVariableIndex('Background_Shape_Error')
			bkg_shape_constraint_val = np.random.randn()*bkg_shape_err
			likelihood.SetGaussianConstraintAbsolute(likelihood.model.variable_list[bkg_shape_idx]['Name'],\
								bkg_shape_constraint_val,\
								bkg_shape_err)


		print('\n\nRunning dataset {}...\n'.format(j))
		likelihood.PrintVariableList()

		print('\nConstraints:')
		for constraint in likelihood.model.constraints:
			print('\t{}'.format(constraint))
		print('\n')

		# Fix the signal parameter to different values, and calculate lambda
		for i in (range(0,num_hypotheses)):

			# Fix the 0nu parameter to a specific hypothesis
			signal_hypothesis = float(i)*2.+0.000001
			signal_idx = likelihood.GetVariableIndex('Bb0n')
			initial_values = np.copy(initial_guess)
			initial_values[signal_idx] = signal_hypothesis
			xvals[i] = signal_hypothesis

				###########################################################################
			# All the exciting stuff happens here!
			lambda_fit_result = likelihood.ComputeLambdaForPositiveSignal(\
									initial_values=initial_values,\
									signal_name='Bb0n',\
									signal_expectation=0.,\
									print_level=1,\
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
			xfit, yfit, crossing, crossing_idx = FindIntersectionByQuadraticInterpolationWilks(\
									xvals[fixed_fit_covar],\
									lambdas[fixed_fit_covar] )

			if DEBUG_PLOTTING:
				plt.clf()
				plt.plot(xfit,np.ones(len(xfit))*2.706,'-b',label='90%CL spline')
				plt.plot(xfit,yfit,'-r',label='Quadratic approx')
				plt.plot( xvals[fixed_fit_converged],\
						lambdas[fixed_fit_converged],\
						'-o',color='tab:blue',label='Actual fits')
				#plt.plot(xvals,lambdas,label='Actual fits')
				plt.plot(crossing,yfit[crossing_idx],'ok')
				plt.xlim(0.,30.)
				plt.ylim(0.,7.)
				plt.xlabel('Num signal')
				plt.ylabel('Lambda')
				plt.legend(loc='upper right')
				plt.savefig('{}/example_fit_curve_DEBUGPLOT_{}.png'.format(output_dir,j),\
									dpi=200,bbox_inches='tight')

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
		#output_row['dataset'] = likelihood.dataset

		output_df_list.append(output_row)


		print('\nDataset {} finished at {:4.4}s'.format(j,time.time()-last_time))
		print('Total time elapsed: {:4.4}s ({:4.4} min)\n\n\n'.format( (time.time() - start_time),\
									(time.time() - start_time)/60. ))
		last_time = time.time()



	##########################################################################
	# Write the output dataframe to a file

	output_df = pd.DataFrame(output_df_list)
	#print(output_df.head())
	print('Saving file to output directory: {}'.format(output_dir))
	if not os.path.exists(output_dir + '/h5/' + date):
		os.makedirs(output_dir + '/h5/' + date)	
	output_df.to_hdf('{}/h5/{}/sens_output_file_90CL_{}_resolution_{}.h5'.format(output_dir,date,iteration_num,resolution),
					 key='df')

	print('Elapsed: {:4.4}s'.format(time.time()-start_time))

