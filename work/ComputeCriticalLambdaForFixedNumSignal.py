# Import sys, then tell python where to find the nEXO-specific classes
import sys
import os
sys.path.append('../modules')


# Check arguments and load inputs
if len(sys.argv) == 4:
	iteration_num = int(sys.argv[1])
	input_num_signal = float(sys.argv[2])
	output_dir = sys.argv[3]
	if not os.path.exists(output_dir):
		sys.exit('\nERROR: path to output_dir does not exist\n')
else:
	print('\n\nERROR: ComputeCriticalLambdaForNumSignal.py requires 3 arguments')
	print('Usage:')
	print('\tpython ConvertExcel2Root_pandas.py <iteration_num> <input_num_signal> </path/to/output/directory/>')
	sys.exit('\n')

# Import useful libraries for analysis
import pandas as pd
import histlite as hl
import numpy as np
from matplotlib import pyplot as plt

# Import the nEXO sensitivity classes
import nEXOFitWorkspace
import nEXOFitModel
import nEXOFitLikelihood

# Import iMinuit
from iminuit import Minuit


# Create the workspace
workspace = nEXOFitWorkspace.nEXOFitWorkspace()
workspace.LoadComponentsTableFromFile('../tables/ComponentsTable_D-005_v25_2020-01-21.h5')
workspace.CreateGroupedPDFs()

# Create the likelihood object
likelihood = nEXOFitLikelihood.nEXOFitLikelihood()
likelihood.AddPDFDataframeToModel(workspace.df_group_pdfs)

initial_values = np.ones(len(likelihood.variable_list))

for i in range(len(likelihood.variable_list)):
	if 'Bb0n' in likelihood.variable_list[i]['Name']:
		# Set signal num_counts to whatever the input value is
		initial_values[i] = input_num_signal
	else:
		# Set background num_counts to the expected values from the groupPDFs.
		initial_values[i] = (likelihood.variable_list[i]['Value'])

likelihood.model_obj.UpdateVariables(initial_values)
likelihood.model_obj.GenerateModelDistribution()
likelihood.AddDataset( likelihood.model_obj.GenerateDataset() )

CONSTRAINTS = True
PAR_LIMITS = True


if PAR_LIMITS:
	# Next, set limits so that none of the PDFs go negative in the fit.
	for var in likelihood.variable_list:
	    if 'Bb0n' in var['Name']:
	        likelihood.SetVariableLimits( var['Name'], \
	                                  lower_limit = 0., \
	                                  upper_limit = 100.)
	    else: 
	        likelihood.SetVariableLimits( var['Name'], \
	                                  lower_limit = 0., \
	                                  upper_limit = var['Value']*10.)


# Increase the step size for the Bb0n variable
likelihood.SetFractionalMinuitError('Num_FullLXeBb0n', 0.01/0.0001)



def NegLogLikelihood(parameter_values):
    return likelihood.ComputeNegLogLikelihood(parameter_values)


# Create Minimizer.
m = Minuit.from_array_func( NegLogLikelihood, \
				initial_values, \
				error = likelihood.GetMinuitErrorTuple(), \
				fix   = likelihood.GetVariableFixTuple(), \
				name  = likelihood.GetVariableNamesTuple(), \
				limit = likelihood.GetVariableLimitsTuple(), \
				errordef = 0.5,\
				print_level = 1)
print('Getting param states:')
try:
	print(m.get_param_states())
except ValueError as e:
	print('ValueError: {}'.format(e))


num_datasets = 5000

output_row = dict()
output_df_list = []

import time

start = time.time()
last_time = start

for j in range(0,num_datasets):
	#print('Running dataset {}'.format(j))
	best_fit_converged = True
	fixedSig_fit_converged = True
	this_lambda = -1.
	output_row = dict()
	likelihood.model_obj.UpdateVariables(initial_values)
	likelihood.model_obj.GenerateModelDistribution()
	likelihood.AddDataset( likelihood.model_obj.GenerateDataset() )

	likelihood.SetAllVariablesFloating()
	likelihood.SetVariableFixStatus('Num_FullTPC_Co60',True)
	
	if CONSTRAINTS:
		rn222_idx = likelihood.GetVariableIndex('Rn222')
		# Fluctuate Rn222 constraint
		rn222_constraint_val = (np.random.randn()*0.1 + 1)*initial_values[rn222_idx]
		# Set Rn222 constraint
		likelihood.SetGaussianConstraintAbsolute(likelihood.variable_list[rn222_idx]['Name'],\
							 rn222_constraint_val, \
	                	                         0.1 * initial_values[rn222_idx])

	print('\n\nRunning dataset {}....\n'.format(j))
	likelihood.PrintVariableList()

	print('\nConstraints:')
	for constraint in likelihood.constraints:
		print('\t{}'.format(constraint))
	print('\n')

	print('\nBest fit:\n')
	m = Minuit.from_array_func( NegLogLikelihood, \
					initial_values, \
					error = likelihood.GetMinuitErrorTuple(), \
					fix   = likelihood.GetVariableFixTuple(), \
					name  = likelihood.GetVariableNamesTuple(), \
					limit = likelihood.GetVariableLimitsTuple(), \
					errordef = 0.5 )
	m.migrad()
	print('Param states after best fit:')
	try:
		print(m.get_param_states())
	except ValueError as e:
		print('ValueError: {}'.format(e))
	nll_best = m.fval
	if not m.get_fmin()['is_valid']:
		best_fit_converged = False
		print('Best fit did not converge.')
	else:
		print('Best fit converged!')
 

	print('\n\nFit with signal value fixed at {:3.3} cts:\n'.format(input_num_signal))

	likelihood.SetVariableFixStatus('Num_FullLXeBb0n',True)

	m = Minuit.from_array_func( NegLogLikelihood, \
			initial_values, \
			error = likelihood.GetMinuitErrorTuple(), \
			fix   = likelihood.GetVariableFixTuple(), \
			name  = likelihood.GetVariableNamesTuple(), \
			limit = likelihood.GetVariableLimitsTuple(), \
			errordef = 0.5 )
	m.migrad()
	print('Param states after fixed-signal fit:')
	try:
		print(m.get_param_states())
	except ValueError as e:
		print('ValueError: {}')
	this_lambda = 2*(m.fval - nll_best)
	if not m.get_fmin()['is_valid']:
		fixedSig_fit_converged = False
		print('Fixed fit did not converge.')
	else:
		print('Fixed fit converged!')

	output_row['num_signal'] = input_num_signal
	output_row['lambda'] = this_lambda
	output_row['best_fit_converged'] = best_fit_converged
	output_row['fixedSig_fit_converged'] = fixedSig_fit_converged
	#output_row['dataset'] = likelihood.dataset

	output_df_list.append(output_row)	
	
	print('Dataset {} finished at {:4.4}s'.format(j,time.time()-last_time))
	last_time = time.time()

output_df = pd.DataFrame(output_df_list)
#print(output_df.head())
print('Saving file to output directory: {}'.format(output_dir))
output_df.to_hdf('{}/critical_lambda_calculation_num_sig_{:3.3}_file_{}.h5'.format(output_dir,input_num_signal,iteration_num),key='df')

print('Elapsed: {:4.4}s'.format(time.time()-start))
