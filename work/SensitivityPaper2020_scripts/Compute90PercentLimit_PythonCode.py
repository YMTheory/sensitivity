# Import sys, then tell python where to find the nEXO-specific classes
import sys
sys.path.append('../modules')
######################################################################
# Load the arguments:
if len(sys.argv) != 4:
	print('\nERROR: incorrect number of arguments.\n')
	print('Usage:')
	print('\tpython Compute90PercentLimit_PythonCode.py ' + \
		'<iteration_num> <critical_lambda_file> <output_directory>\n\n')
	sys.exit()

iteration_num = int(sys.argv[1])
critical_lambda_file = sys.argv[2]
output_dir = sys.argv[3]
#####################################################################

##########################################################################
##########################################################################
def FindIntersectionByQuadraticInterpolation( xvals, yvals, SplineFunction ):

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
		yspline = SplineFunction(xfit)
		crossing_idx = np.where( (yfit - yspline)>0. )[0][0]
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
from scipy.interpolate import LSQUnivariateSpline

# Import the nEXO sensitivity classes
import nEXOFitWorkspace
import nEXOFitModel
import nEXOFitLikelihood


# Get the critical lambda data points, fit them to a spline curve.
critical_lambda_data = np.genfromtxt(critical_lambda_file)
spline_xn = np.array([1., 5., 7., 10., 20., 30, 48.5]) # defines the locations of the knots
SplineFunc = LSQUnivariateSpline(critical_lambda_data[:,0],critical_lambda_data[:,1],t = spline_xn,k=3)

# Create the workspace
workspace = nEXOFitWorkspace.nEXOFitWorkspace()
workspace.LoadComponentsTableFromFile('../tables/ComponentsTable_D-005_v25_2020-01-21.h5')
workspace.CreateGroupedPDFs()


# Create the likelihood object
likelihood = nEXOFitLikelihood.nEXOFitLikelihood()
likelihood.AddPDFDataframeToModel(workspace.df_group_pdfs)
initial_guess = np.ones(len(likelihood.variable_list))

for i in range(len(likelihood.variable_list)):
	initial_guess[i] = (likelihood.variable_list[i]['Value'])

likelihood.model_obj.UpdateVariables(initial_guess)
likelihood.model_obj.GenerateModelDistribution()
likelihood.AddDataset( likelihood.model_obj.GenerateDataset() )




PAR_LIMITS = True
DEBUG_PLOTTING = False

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

likelihood.SetFractionalMinuitInputError('Num_FullLXeBb0n', 0.01/0.0001)


##########################################################################
# Here's where the calculation loop begins.

num_datasets = 100
num_hypotheses = 20
xvals = np.array([]) #np.linspace(0.,40.*(1.-num_hypotheses),num_hypotheses)

lambdas = np.zeros(num_hypotheses)
num_iterations = np.zeros(num_hypotheses)
best_fit_converged = True
#fixed_fit_converged = np.ones(num_hypotheses,dtype=bool)
#fixed_fit_covar = np.ones(num_hypotheses, dtype=bool)
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

	likelihood.model_obj.UpdateVariables(initial_guess)
	likelihood.model_obj.GenerateModelDistribution()
	likelihood.AddDataset( likelihood.model_obj.GenerateDataset() )
	
	likelihood.SetAllVariablesFloating()

        # Fix the Co60 parameters
	likelihood.SetVariableFixStatus('Num_FullTPC_Co60',True)

	RN_CONSTRAINTS=True
	if RN_CONSTRAINTS:
		rn222_idx = likelihood.GetVariableIndex('Rn222')
		# Fluctuate Rn222 constraint
		rn222_constraint_val = (np.random.randn()*0.1 + 1)*initial_guess[rn222_idx]
		# Set Rn222 constraint
		likelihood.SetGaussianConstraintAbsolute(likelihood.variable_list[rn222_idx]['Name'],\
							 rn222_constraint_val, \
	                	                         0.1 * initial_guess[rn222_idx])

	print('\n\nRunning dataset {}...\n'.format(j))
	likelihood.PrintVariableList()	

	print('\nConstraints:')
	for constraint in likelihood.constraints:
		print('\t{}'.format(constraint))
	print('\n')


	# Find the global best fit first
	print('\n\nBefore fit starts:')
	likelihood.PrintVariableList()
	best_fit_iterations = likelihood.CreateAndRunMinuitFit( initial_guess, print_level=1 )
	print('After fit ends:')
	likelihood.PrintVariableList()
	nll_best = likelihood.fitter.fval
	if not likelihood.fitter.get_fmin()['is_valid']:
		best_fit_converged = False
		print('Best fit did not converge.')
	else:
		print('Best fit converged!')
	

	# Next, fix the signal parameter to specific hypotheses
	for i in (range(0,num_hypotheses)):

		# Fix the 0nu parameter to a specific hypothesis
		signal_idx = likelihood.GetVariableIndex('Bb0n')
		initial_values = np.copy(initial_guess)
		initial_values[signal_idx] = float(i)*2.+0.000001
		xvals[i] = float(i)*2.+0.000001
		likelihood.SetVariableFixStatus('Num_FullLXeBb0n',True)    
		print('\n\nHypothesis {} ({:3.3} signal counts)'.format(i,initial_values[signal_idx]))
		#print('\n\nBefore fit starts:')
		#likelihood.PrintVariableList()

		num_iterations[i] = likelihood.CreateAndRunMinuitFit( initial_values, print_level=1 )

		lambdas[i] = 2*(likelihood.fitter.fval - nll_best)

		if not likelihood.fitter.get_fmin()['is_valid']:
			fixed_fit_converged = np.append(fixed_fit_converged,np.array([False]))
		else:
			fixed_fit_converged = np.append(fixed_fit_converged,np.array([True]))
		if not likelihood.fitter.get_fmin()['has_accurate_covar']:
			fixed_fit_covar = np.append(fixed_fit_covar,np.array([False]))
		else:
			fixed_fit_covar = np.append(fixed_fit_covar,np.array([True]))


		print('After fit ends:')
		#print(likelihood.fitter.get_fmin())
			#print('{}'.format(key))
		likelihood.PrintVariableList()
		print('\n')

	# Next, find the hypothesis value for which lambda crosses the critical lambda curve 
	if best_fit_converged:
		xfit, yfit, crossing, crossing_idx = FindIntersectionByQuadraticInterpolation(\
								xvals[fixed_fit_converged],\
								lambdas[fixed_fit_converged],\
								SplineFunc )
	
		if DEBUG_PLOTTING:
			plt.clf()
			plt.plot(xfit,SplineFunc(xfit),'-b',label='90%CL spline')
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
	output_row['best_fit_converged'] = best_fit_converged
	output_row['fixed_fit_converged'] = fixed_fit_converged
	output_row['fixed_fit_acc_covar'] = fixed_fit_covar
	output_row['90CL_crossing'] = crossing
	output_row['num_iterations'] = num_iterations
	output_row['best_fit_iterations'] = best_fit_iterations
	#output_row['dataset'] = likelihood.dataset

	output_df_list.append(output_row)
	
	
	print('\nDataset {} finished at {:4.4}s.'.format(j,time.time()-last_time))
	print('Total time elapsed: {:4.4}s ({:4.4} min)\n\n\n'.format( (time.time() - start_time),\
								(time.time() - start_time)/60. ))
	last_time = time.time()



##########################################################################
# Write the output dataframe to a file

output_df = pd.DataFrame(output_df_list)
#print(output_df.head())
print('Saving file to output directory: {}'.format(output_dir))
output_df.to_hdf('{}/sens_output_file_90CL_{}_00.h5'.format(output_dir,iteration_num),key='df')

print('Elapsed: {:4.4}s'.format(time.time()-start_time))





#	mask = ((yvals>0.5)&(yvals<6.)&fixed_fit_converged)[1:]
#	x = np.linspace(0.,40.,800)
#	if np.sum(mask) == 0 or len( (xvals[1:])[mask] ) == 0:
#			continue
#
#	# If the crossing is low, include the first bin; otherwise
#	# the fit seems to be biased and get the wrong crossing point.
#	if crossing < 5. and yfit[0] < yspline[0]:
#		mask = ((yvals>0.)&(yvals<6.)&fixed_fit_converged)
#		if len( (xvals[mask]) ) == 0:
#			crossing = -1. 
#			continue
#		else:
#			try:
#				p = np.polyfit( (xvals)[mask], (yvals)[mask], 2.)
#				yfit = p[0]*x**2 + p[1]*x + p[2]
#			except np.RankWarning:
#				p = np.polyfit( (xvals)[mask], (yvals)[mask], 1.)
#				yfit = p[0]*x + p[1]
#				#yfit = p[0]*x**2 + p[1]*x + p[2]
#			crossing_idx = np.where( (yfit - yspline)>0. )[0][0]
#			crossing = x[crossing_idx]
#	# If the fit goes ABOVE the spline at 0, exclude the first THREE bins
#	# to ensure we're looking at only the upper limit crossing. Otherwise, 
#	# the "crossing" will be reconstructed at 0 given the above code
#	elif crossing < 2. and yfit[0] > yspline[0]:
#		mask = ((yvals>1.)&(yvals<6.)&fixed_fit_converged)[3:]
#		if len( (xvals[3:])[mask] ) > 0:
#			crossing = -1.
#			continue
#		else:
#			try:
#				p = np.polyfit( (xvals[3:])[mask], (yvals[3:])[mask], 2.)
#				yfit = p[0]*x**2 + p[1]*x + p[2]
#			except np.RankWarning:
#				p = np.polyfit( (xvals[3:])[mask], (yvals[3:])[mask], 1.)
#			yfit = p[0]*x**2 + p[1]*x + p[2]
#			crossing_idx = np.where( (yfit - yspline)>0. )[0][0]
#			crossing = x[crossing_idx]


