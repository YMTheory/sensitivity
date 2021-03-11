# Import sys, then tell python where to find the nEXO-specific classes
import sys
sys.path.append('../../../modules')
######################################################################
# Load the arguments:
if len(sys.argv) != 10:
	print('\nERROR: incorrect number of arguments.\n')
	print('Usage:')
	print('\tpython Compute90PercentLimit_PythonCode.py ' + \
		'<job_id_num> <bkg_shape_err_parameter> ' + \
                '<num_datasets_to_generate> <input_table> <output_directory> ' + \
                '<bb0n_count> <livetime>\n\n')
	sys.exit()

job_id_num = int(sys.argv[1])
bkg_shape_err = float(sys.argv[2])
num_datasets = int(sys.argv[3])
input_table = sys.argv[4]
output_dir = sys.argv[5]
bb0n_count = int(sys.argv[6])
livetime = float(sys.argv[7])
yaml_card = sys.argv[8]
xe137_scale_factor = float(sys.argv[9])


num_hypotheses = 1
step_hypothesis = 2

#####################################################################
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
INCLUDE_BACKGROUND_SHAPE_ERROR = False
PAR_LIMITS = True

# Create the workspace
workspace = nEXOFitWorkspace.nEXOFitWorkspace(yaml_card)
workspace.LoadComponentsTableFromFile(input_table)

# Scale the Xe137 component. Note that this is done at the ComponentsTable level,
# since (unlike Rn222) the Xe137 doesn't have its own group.
for index, row in workspace.df_components.iterrows():
	if 'Xe137' in row['PDFName']:
		workspace.df_components.loc[index,'SpecActiv'] = xe137_scale_factor * row['SpecActiv']
		workspace.df_components.loc[index,'SpecActivErr'] = xe137_scale_factor * row['SpecActivErr']

workspace.livetime = livetime * 365. * 24. * 60. * 60.
workspace.CreateGroupedPDFs()

# Create the likelihood object
likelihood = nEXOFitLikelihood.nEXOFitLikelihood()
likelihood.AddPDFDataframeToModel(workspace.df_group_pdfs, workspace.histogram_axis_names)


if INCLUDE_BACKGROUND_SHAPE_ERROR:
	likelihood.model.IncludeBackgroundShapeVariableInFit(True)

initial_guess = likelihood.GetVariableValues()

# Get the initial values;
initial_values = likelihood.GetVariableValues()

# Set the BB0n num_signal to the user-provided input
initial_values[likelihood.GetVariableIndex('Bb0n')] = bb0n_count

# Update the model in the likelihood object
likelihood.model.UpdateVariables(initial_guess)
likelihood.model.GenerateModelDistribution()
 
# Add an initial dataset
likelihood.AddDataset(likelihood.model.GenerateDataset())

 
# Set limits so that none of the PDFs can go negative in the fit.
# (except the signal PDF)
if PAR_LIMITS:
	for var in likelihood.model.variable_list:
		if 'Bb0n' in var['Name']:
			likelihood.SetVariableLimits( var['Name'], lower_limit = -15.0, upper_limit = 100.)
		elif 'Background_Shape_Error' in var['Name']:
			likelihood.SetVariableLimits( var['Name'], lower_limit = -100.0, upper_limit = 100.)
		elif 'Co60' in var['Name']:
			likelihood.SetVariableLimits( var['Name'], lower_limit = 0.0, upper_limit = var['Value']*10.)
		else:
			likelihood.SetVariableLimits( var['Name'], lower_limit = var['Value']*0.1, upper_limit = var['Value']*10.0)

likelihood.SetFractionalMinuitInputError('Num_FullLXeBb0n', 0.01/0.0001)


##########################################################################
# Here's where the calculation loop begins.
workspace.SetHandlingOfRadioassayData(fluctuate=True)

output_df_list = []

import time

start_time = time.time()
last_time = start_time

for j in range(0,num_datasets):
	# Initialize (or reset) all my output variables.

	num_iterations = np.ones(num_hypotheses)
	lambdas = np.zeros(num_hypotheses)
	xvals = np.zeros(num_hypotheses)
	fixed_fit_converged = np.array([],dtype=bool)
	fixed_fit_covar = np.array([],dtype=bool)
	fixed_fit_parameters = np.array([],dtype=bool)
	fixed_fit_errors = np.array([],dtype=bool)
	output_row = dict()

	workspace.CreateGroupedPDFs()
	likelihood.AddPDFDataframeToModel(workspace.df_group_pdfs, workspace.histogram_axis_names, replace_existing_variables=False)

	# likelihood.model.UpdateVariables(initial_guess)
	sig_idx = likelihood.model.GetVariableIndexByName('FullLXeBb0n')
	likelihood.model.variable_list[sig_idx]['Value'] = bb0n_count

	if INCLUDE_BACKGROUND_SHAPE_ERROR: 
		shape_idx = likelihood.model.GetVariableIndexByName('Background_Shape_Error')
		likelihood.model.variable_list[shape_idx]['Value'] = 0

	initial_guess = likelihood.GetVariableValues()
	likelihood.model.GenerateModelDistribution()
	likelihood.AddDataset(likelihood.model.GenerateDataset())
	
	
	likelihood.SetAllVariablesFloating()

    #Fix the Co60 parameter
	# likelihood.SetVariableFixStatus('Num_FullTPC_Co60', True)


	RN_CONSTRAINTS=True
	if RN_CONSTRAINTS:
		xe137_idx = likelihood.GetVariableIndex('Rn222')
		# Fluctuate Rn222 constraint
		rn222_constraint_val = (np.random.randn() * 0.1 + 1) * initial_guess[xe137_idx]
		# Set Rn222 constraint
		likelihood.SetGaussianConstraintAbsolute(likelihood.model.variable_list[xe137_idx]['Name'], rn222_constraint_val, 0.1 * initial_guess[xe137_idx])

	B8_CONSTRAINTS=True
	if B8_CONSTRAINTS:
		b8_index = likelihood.GetVariableIndex('B8')
		# Fluctuate B8nu constraint
		b8_constraint_val = (np.random.randn() * 0.1 + 1) * initial_guess[b8_index]
		# Set B8nu constraint
		likelihood.SetGaussianConstraintAbsolute(likelihood.model.variable_list[b8_index]['Name'], b8_constraint_val, 0.1 * initial_guess[b8_index])

			
	if INCLUDE_BACKGROUND_SHAPE_ERROR:
		bkg_shape_idx = likelihood.GetVariableIndex('Background_Shape_Error')
		bkg_shape_constraint_val = np.random.randn()*bkg_shape_err
		likelihood.SetGaussianConstraintAbsolute(likelihood.model.variable_list[bkg_shape_idx]['Name'], bkg_shape_constraint_val, bkg_shape_err)


	print('\n\nRunning dataset {}...\n'.format(j))
	likelihood.PrintVariableList()	

	print('\nConstraints:')
	for constraint in likelihood.model.constraints:
		print('\t{}'.format(constraint))
	print('\n')

	# Fix the signal parameter to different values, and calculate lambda
	for i in (range(0,num_hypotheses)):

		# Fix the 0nu parameter to a specific hypothesis
		signal_hypothesis = float(i)*step_hypothesis + 0.000001
		signal_idx = likelihood.GetVariableIndex('Bb0n')
		initial_values = np.copy(initial_guess)
		initial_values[signal_idx] = signal_hypothesis
		xvals[i] = signal_hypothesis
		
		print(initial_values)

	    ###########################################################################
		# All the exciting stuff happens here!
		lambda_fit_result = likelihood.ComputeLambdaForPositiveSignal(\
								initial_values=initial_values,\
								signal_name='Bb0n',\
								signal_expectation=0.0,\
								print_level=1,\
								fixed_fit_signal_value=signal_hypothesis)
	    ###########################################################################
								
		print(lambda_fit_result['best_fit_parameters'])			

		# Store all the important quantities	
		lambdas[i] = lambda_fit_result['lambda']
		fixed_fit_converged = np.append(fixed_fit_converged, lambda_fit_result['fixed_fit_converged'])
		fixed_fit_covar = np.append(fixed_fit_covar, lambda_fit_result['fixed_fit_covar'])
		fixed_fit_parameters = np.append(fixed_fit_parameters, lambda_fit_result['fixed_fit_parameters'])
		fixed_fit_errors = np.append(fixed_fit_errors, lambda_fit_result['fixed_fit_errors'])
		num_iterations[i] = lambda_fit_result['fixed_fit_iterations']

		print('After fit ends:')
		likelihood.PrintVariableList()
		print('\n')


	output_row['num_signal'] = xvals
	output_row['lambda'] = lambdas
	output_row['fixed_fit_converged'] = fixed_fit_converged
	output_row['fixed_fit_acc_covar'] = fixed_fit_covar
	output_row['fixed_fit_parameters'] = fixed_fit_parameters
	output_row['fixed_fit_errors'] = fixed_fit_errors
	output_row['num_iterations'] = num_iterations
	output_row['best_fit_converged'] = lambda_fit_result['best_fit_converged']
	output_row['best_fit_covar'] = lambda_fit_result['best_fit_covar']
	output_row['best_fit_iterations'] = lambda_fit_result['best_fit_iterations']
	output_row['best_fit_parameters'] = lambda_fit_result['best_fit_parameters']
	output_row['best_fit_errors'] = lambda_fit_result['best_fit_errors']
	output_row['input_parameters'] = initial_guess
	output_df_list.append(output_row)
	
	print('\nDataset {} finished at {:4.4}s'.format(j,time.time()-last_time))
	print('Total time elapsed: {:4.4}s ({:4.4} min)\n\n\n'.format( (time.time() - start_time), (time.time() - start_time)/60. ))
	last_time = time.time()



##########################################################################
# Write the output dataframe to a file

output_df = pd.DataFrame(output_df_list)
#print(output_df.head())
print('Saving file to output directory: {}'.format(output_dir))
output_df.to_hdf('{}/DiscoveryPotential_bb0n_{}ct_{}yrs_g{:03}_{:03}.h5'.format(output_dir, bb0n_count, livetime, xe137_scale_factor, job_id_num ),key='df')

print('Elapsed: {:4.4}s'.format(time.time()-start_time)) 



