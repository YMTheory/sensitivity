# Import sys, then tell python where to find the nEXO-specific classes
import sys
import os
sys.path.append('../../modules')


######################################################################
# Check arguments and load inputs
if len(sys.argv) == 6:
	iteration_num = int(sys.argv[1])
	input_num_signal = float(sys.argv[2])
	num_datasets = int(sys.argv[3])
	input_table = sys.argv[4]
	output_dir = sys.argv[5]
	if not os.path.exists(output_dir):
		sys.exit('\nERROR: path to output_dir does not exist\n')
else:
	print('\n\nERROR: ComputeCriticalLambdaForNumSignal.py requires 7 arguments')
	print('Usage:')
	print('\tpython ComputeCriticalLambdaForNumSignal.py ' + \
		'<iteration_num> <input_num_signal> <num_datasets_to_generate> <input_table> </path/to/output/directory/> ')
	sys.exit('\n')
######################################################################

mass_str = (input_table.split('_')[-1]).split('k')[0][5:]
print('\nFiducial mass: {} kg\n'.format(mass_str))
mass_int = int(mass_str)


# Import useful libraries for analysis
import pandas as pd
import histlite as hl
import numpy as np
from matplotlib import pyplot as plt
import copy

# Import the nEXO sensitivity classes
import nEXOFitWorkspace
import nEXOFitModel
import nEXOFitLikelihood

# Import iMinuit
from iminuit import Minuit


# Create the workspace
workspace = nEXOFitWorkspace.nEXOFitWorkspace(config='/p/vast1/nexo/sensitivity2020/pdfs/'+\
                   'config_files/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1.yaml')
workspace.LoadComponentsTableFromFile(input_table)
workspace.CreateGroupedPDFs()


# Create the likelihood object
likelihood = nEXOFitLikelihood.nEXOFitLikelihood()
likelihood.AddPDFDataframeToModel(workspace.df_group_pdfs, workspace.histogram_axis_names)


# Get the initial values; set the BB0n num_signal to the user-provided input
sig_idx = likelihood.model.GetVariableIndexByName('Bb0n')
likelihood.model.variable_list[ sig_idx ]['Value'] = input_num_signal


initial_values = likelihood.GetVariableValues()

# Initialize the model
likelihood.model.UpdateVariables(initial_values)
likelihood.model.GenerateModelDistribution()

likelihood.AddDataset( likelihood.model.GenerateDataset() )


CONSTRAINTS = True
PAR_LIMITS = True

if PAR_LIMITS:
	# Next, set limits so that none of the PDFs go negative in the fit.
	for var in likelihood.model.variable_list:
	    if 'Bb0n' in var['Name']:
	        likelihood.SetVariableLimits( var['Name'], \
	                                  lower_limit = 0., \
	                                  upper_limit = 100.)
	    elif 'Co60' in var['Name']:
                likelihood.SetVariableLimits( var['Name'], \
                                          lower_limit = 0., \
                                          upper_limit = var['Value']*10.)
	    else: 
	        likelihood.SetVariableLimits( var['Name'], \
	                                  lower_limit = 0.1*var['Value'], \
	                                  upper_limit = var['Value']*10.)


# Increase the step size for the Bb0n variable
likelihood.SetFractionalMinuitInputError('Num_FullLXeBb0n', 0.01/0.0001)

workspace.SetHandlingOfRadioassayData( fluctuate=True )

output_row = dict()
output_df_list = []

import time

start_time = time.time()
last_time = start_time

for j in range(0,num_datasets):

	best_fit_converged = True
	fixedSig_fit_converged = True
	this_lambda = -1.
	output_row = dict()
	best_fit_parameters = None
	best_fit_errors = None
	fixed_fit_parameters = None
	fixed_fit_errors = None


	#likelihood.model.UpdateVariables(initial_values)


	# Redo the grouping, which fluctuates the radioassay values within their uncertainties.
	workspace.CreateGroupedPDFs()
	likelihood.AddPDFDataframeToModel( workspace.df_group_pdfs, \
                                           workspace.histogram_axis_names, \
                                           replace_existing_variables=False )

	sig_idx = likelihood.model.GetVariableIndexByName('Bb0n')
	likelihood.model.variable_list[ sig_idx ]['Value'] = input_num_signal

	# Need to regenerate the distribution since we've changed Num_Bb0n
	likelihood.model.GenerateModelDistribution()
	likelihood.AddDataset( likelihood.model.GenerateDataset() )

	# Reset the num signal, also save input values as dict
	input_parameters = dict()
	for var in likelihood.model.variable_list:
		if 'Bb0n' in var['Name']:
			input_parameters[ var['Name'] ] = input_num_signal
		else:
			input_parameters[ var['Name'] ] = float(var['Value'])

	likelihood.SetAllVariablesFloating()

	#likelihood.SetVariableFixStatus('Num_FullTPC_Co60',True)
	
	if CONSTRAINTS:
		rn222_idx = likelihood.GetVariableIndex('Rn222')
		# Fluctuate Rn222 constraint
		rn222_constraint_val = (np.random.randn()*0.1 + 1)*initial_values[rn222_idx]
		# Set Rn222 constraint
		likelihood.SetGaussianConstraintAbsolute(likelihood.model.variable_list[rn222_idx]['Name'],\
							 rn222_constraint_val, \
	                	                         0.1 * initial_values[rn222_idx])
		b8_idx = likelihood.GetVariableIndex('B8')
		# Fluctuate B8nu constraint
		b8_constraint_val = (np.random.randn()*0.1 + 1)*initial_values[b8_idx]
		# Set B8nu constraint
		likelihood.SetGaussianConstraintAbsolute(likelihood.model.variable_list[b8_idx]['Name'],\
							 b8_constraint_val, \
	                	                         0.1 * initial_values[b8_idx])


	print('\n\nRunning dataset {}....\n'.format(j))
	likelihood.PrintVariableList()

	print('\nConstraints:')
	for constraint in likelihood.model.constraints:
		print('\t{}'.format(constraint))
	print('\n')
	
	###########################################################################
	# All the exciting stuff happens here!
	lambda_fit_result = likelihood.ComputeLambdaForPositiveSignal(\
							initial_values=initial_values,\
							signal_name='Bb0n',\
							signal_expectation=input_num_signal,\
							print_level=1 )
	###########################################################################

	output_row['num_signal']           = input_num_signal
	output_row['lambda']               = lambda_fit_result['lambda']
	output_row['best_fit_converged']   = lambda_fit_result['best_fit_converged']
	output_row['best_fit_covar']       = lambda_fit_result['best_fit_covar']
	output_row['best_fit_iterations']  = lambda_fit_result['best_fit_iterations']
	output_row['fixed_fit_converged']  = lambda_fit_result['fixed_fit_converged']
	output_row['fixed_fit_covar']      = lambda_fit_result['fixed_fit_covar']
	output_row['fixed_fit_iterations'] = lambda_fit_result['fixed_fit_iterations']
	output_row['best_fit_parameters']  = lambda_fit_result['best_fit_parameters']
	output_row['best_fit_errors']      = lambda_fit_result['best_fit_errors']
	output_row['fixed_fit_parameters'] = lambda_fit_result['fixed_fit_parameters']
	output_row['fixed_fit_errors']     = lambda_fit_result['fixed_fit_errors']
	output_row['input_parameters']     = input_parameters
	output_row['best_fit_nll']         = lambda_fit_result['best_fit_nll']
	#output_row['dataset'] = likelihood.dataset

	output_df_list.append(output_row)	
	print('Variable values at end of loop:')
	likelihood.PrintVariableList()
	
	print('Dataset {} finished at {:4.4}s'.format(j,time.time()-last_time))
	last_time = time.time()

output_df = pd.DataFrame(output_df_list)
#print(output_df.head())
print('Saving file to output directory: {}'.format(output_dir))
output_df.to_hdf('{}/critical_lambda_calculation_fidvolstudy_{}_D-024_num_sig_{:06.6}_file_{}.h5'.format(\
			output_dir, mass_str, input_num_signal, iteration_num),key='df')

print('Elapsed: {:4.4}s'.format(time.time()-start_time))