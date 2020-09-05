# Import sys, then tell python where to find the nEXO-specific classes
import sys
sys.path.append('../../modules')

######################################################################
# Load the arguments:
if len(sys.argv) != 7:
	print('\nERROR: incorrect number of arguments.\n')
	print('Usage:')
	print('\tpython Compute90PercentLimit_PythonCode.py ' + \
		'<job_id_num> <bkg_shape_err_parameter> ' + \
                '<num_datasets_to_generate> <input_table> <output_directory> ' + \
                '<xe137_scale_factor>\n\n')
	sys.exit()

job_id_num = int(sys.argv[1])
bkg_shape_err = float(sys.argv[2])
num_datasets = int(sys.argv[3])
input_table = sys.argv[4]
output_dir = sys.argv[5]
xe137_scale_factor = float(sys.argv[6])

#####################################################################

##########################################################################
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

	xfit = np.linspace(0.,60.,1200)

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
INCLUDE_BACKGROUND_SHAPE_ERROR = True
PAR_LIMITS = True
DEBUG_PLOTTING = False


# Create the workspace
workspace = nEXOFitWorkspace.nEXOFitWorkspace('./config/TUTORIAL_config.yaml')
workspace.LoadComponentsTableFromFile( input_table )
workspace.CreateGroupedPDFs()

# Define the ROI within the workspace
roi_dict = { 'SS/MS':         [0.,1.],
             'Energy (keV)':  [ 2428.89, 2486.77 ], 
             'Standoff (mm)': [ 120., 650. ] }
workspace.DefineROI( roi_dict )



# Create the likelihood object
likelihood = nEXOFitLikelihood.nEXOFitLikelihood()
likelihood.AddPDFDataframeToModel(workspace.df_group_pdfs)


if INCLUDE_EFFICIENCY_ERROR:
	likelihood.model.IncludeSignalEfficiencyVariableInFit(True)
if INCLUDE_BACKGROUND_SHAPE_ERROR:
	likelihood.model.IncludeBackgroundShapeVariableInFit(True)


initial_guess = likelihood.GetVariableValues()

# Scale the Xe13 component according to the input value
xe137_idx = likelihood.model.GetVariableIndexByName('Xe137')
initial_guess[xe137_idx] *= xe137_scale_factor

# Update the model in the likelihood object
likelihood.model.UpdateVariables(initial_guess)
likelihood.model.GenerateModelDistribution()

# Print out the number of events in the ROI
total_bkg_in_roi = likelihood.model.GetIntegralInBinRange( workspace.GetROIBinIndices() )
print('\n****************************************************************************************')
print('Variable list after scaling Xe137 component:')
likelihood.PrintVariableList()
print('\n\nTotal ROI background: {:4.4} events\n'.format( \
      likelihood.model.GetIntegralInBinRange( workspace.GetROIBinIndices() ) ))
print('Contribution from each component in ROI:')
for component in likelihood.model.variable_list:
    if 'Shape' in component['Name']: continue
    num_counts_in_roi = likelihood.model.GetComponentIntegralInBinRange( \
                                            component['Name'], workspace.GetROIBinIndices() )
    print('{:<20}\t{:>10.4}\t{:>10.4}%'.format( component['Name']+':', num_counts_in_roi, \
                                              int(1000*num_counts_in_roi/total_bkg_in_roi)/10. ) )
print('****************************************************************************************\n')

# Add an initial dataset
likelihood.AddDataset( likelihood.model.GenerateDataset() )

# Set limits so that none of the PDFs can go negative in the fit.
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

num_hypotheses = 30
xvals = np.array([])

lambdas = np.zeros(num_hypotheses)
num_iterations = np.zeros(num_hypotheses)
best_fit_converged = True
crossing = -1
output_df_list = []

import time

start_time = time.time()
last_time = start_time

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

	likelihood.model.UpdateVariables(initial_guess)
	likelihood.model.GenerateModelDistribution()
	likelihood.AddDataset( likelihood.model.GenerateDataset() )
	
	likelihood.SetAllVariablesFloating()

        # Fix the Co60 parameter
	likelihood.SetVariableFixStatus('Num_FullTPC_Co60',True)

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
	print('Total time elapsed: {:4.4}s ({:4.4} min)\n\n\n'.format( \
                                                                (time.time() - start_time),\
								(time.time() - start_time)/60. ))
	last_time = time.time()



##########################################################################
# Write the output dataframe to a file

output_df = pd.DataFrame(output_df_list)
#print(output_df.head())
print('Saving file to output directory: {}'.format(output_dir))
output_df.to_hdf('{}/sens_output_file_xe137study_{:0>3.3}x_90CL_{:03}.h5'.format(\
                         output_dir, xe137_scale_factor, job_id_num ),\
                         key='df')

print('Elapsed: {:4.4}s'.format(time.time()-start_time))



