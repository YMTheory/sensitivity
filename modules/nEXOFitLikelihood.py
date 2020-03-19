import numpy as np
import histlite as hl
import pandas as pd
import nEXOFitModel
import copy
import sys

import matplotlib.pyplot as plt
from cycler import cycler

from iminuit import Minuit

class nEXOFitLikelihood:

   def __init__( self ):
       # nEXOFitLikelihood is a class that will generate a likelihood function
       # given some arbitrary binned likelihoods and binned PDFs. Both PDFs and
       # datasets are assumed to be in the form of Histlite histograms, and 
       # we assume we're performing a binned maximum likelihood fit.
       self.pdfs = {}
       self.dataset = None
       self.model_distribution = None
       #self.variable_list = []
       self.initial_values = []
       self.model = nEXOFitModel.nEXOFitModel()
       self.nll = np.nan
       self.nll_offset = 0. # -604307041. # This is a legacy number. Offset should be
                                     # recomputed each time you generate a new dataset. 
       self.nll_extra_offset = 0. # This is used to ensure that the minimum is close to zero
                                  # on the second iteration of the fit.
       #self.constraints = []
       self.fitter = None

   ##############################################################################################
   def AddDataset( self, input_dataset ):
       self.dataset = input_dataset
       if self.model_distribution is not None:
          self.nll_offset = 0.
          self.SetInitialOffset() 

   ##############################################################################################
   def AddPDFDataframeToModel( self, df_pdfs ):
       self.model.AddPDFsFromDataframe( df_pdfs )
       self.model_distribution = self.model.GenerateModelDistribution()

   ##############################################################################################
   def ComputeNegLogLikelihood( self, var_values ):
     
       fast = False
       self.model.UpdateVariables( var_values )
       self.model_distribution = self.model.GenerateModelDistribution(fast=fast)
       if not fast:
          self.model_distribution = self.model_distribution.values
       #self.variable_list = self.model.variable_list
      
       # Here we do a binned negative log-likelihood, assuming
       # each bin is an independent, poisson-distributed variable
       mask = (self.model_distribution > 0.)
       self.nll = np.sum( self.model_distribution[mask]  - \
                  self.dataset.values[mask] * np.log( self.model_distribution[mask] ) ) -\
                  self.nll_offset - self.nll_extra_offset 


       # Here we apply the constraints (if any) to the variables.
       for constraint in self.model.constraints:
           self.nll += ( var_values[ constraint['Index'] ] - constraint['Value'] )**2 / \
                       ( np.sqrt(2) * constraint['Width']**2 )

       return self.nll

   ##############################################################################################
   def CreateAndRunMinuitFit( self, initial_values=None, print_level=0 ):

       if initial_values is None:
          initial_values = self.GetVariableValues()             

       self.fitter = Minuit.from_array_func( self.ComputeNegLogLikelihood, \
                                   np.copy(initial_values),\
                                   error = self.GetMinuitInputErrorTuple(), \
                                   fix   = self.GetVariableFixTuple(), \
                                   name  = self.GetVariableNamesTuple(), \
                                   limit = self.GetVariableLimitsTuple(), \
                                   errordef = 0.5,\
                                   print_level = print_level)
       self.fitter.migrad()
       #self.PrintVariableList()
       num_iterations = 1
       while not (self.fitter.get_fmin()['is_valid']\
             and self.fitter.get_fmin()['has_accurate_covar']):    
                 if num_iterations > 9:
                    break
                 if print_level > 0:
                    print('Fit iteration {}'.format(num_iterations+1))

                 #########################################################
                 # Fluctuate the inputs on the second pass, but make sure the
                 # 'fixed' variables stay the same
                 fluctuated_input_values = self.GetVariableValues() + \
                                     np.random.randn(len(initial_values)) * \
                                     np.sqrt( self.GetVariableValues() )
                 for i in range(0,len(self.model.variable_list)):
                          if self.model.variable_list[i]['IsFixed']:
                             fluctuated_input_values[i] = initial_values[i]
                 ##########################################################
          
                 self.fitter = Minuit.from_array_func( self.ComputeNegLogLikelihood, \
                                        np.copy(fluctuated_input_values),\
                                        error = self.GetMinuitInputErrorTuple(), \
                                        fix   = self.GetVariableFixTuple(), \
                                        name  = self.GetVariableNamesTuple(), \
                                        limit = self.GetVariableLimitsTuple(), \
                                        errordef = 0.5,\
                                        print_level = print_level)
                 self.fitter.migrad()
                 num_iterations += 1
                 #self.PrintVariableList()

       fit_errors = self.fitter.errors
       for var in self.model.variable_list:
           var['FitError'] = fit_errors[ var['Name'] ]
           

       return self.fitter.get_fmin()['is_valid'],\
              self.fitter.get_fmin()['has_accurate_covar'],\
              num_iterations 

   ##############################################################################################


   ##############################################################################################
   def GetVariableValues( self ):
       # Returns a numpy array with the variable values
       array = np.zeros(len(self.model.variable_list))
       for i in range(0,len(self.model.variable_list)):
           array[i] = np.copy(self.model.variable_list[i]['Value'])
       return array

   ##############################################################################################
   def SetVariableLimits( self, var_name, upper_limit=None, lower_limit=None):
       var_idx = self.GetVariableIndex( var_name )
       self.model.variable_list[var_idx]['Limits'] = (lower_limit,upper_limit)

   ##############################################################################################
   def PrintVariableList( self ):
       print('{:<21} {:<12} {:<9} {:<10} {:<13} {:<14} {:<13}'.format('Variable name:','Value:',\
                                                        'IsFixed:','FitError','InputError:','IsConstrained:','Limits:'))
       for var in self.model.variable_list:
           print('{:<21} {:<12.4} {:<9} {:<10.4} {:<13.4} {:<14} ({:4},{:4})'.format(var['Name'], \
                                                                      var['Value'],\
                                                                  str(var['IsFixed']),\
                                                                  str(var['FitError']),\
                                                                      var['MinuitInputError'],\
                                                                  str(var['IsConstrained']),\
                                                                  str(var['Limits'][0]),str(var['Limits'][1]) ))

   ##############################################################################################
   def GetVariableIndex( self, var_name ):
       return self.model.GetVariableIndexByName( var_name )

   ##############################################################################################
   def SetInitialOffset( self ):
       if self.model.initial_variable_list == []:
          print('ERROR: Attempting to compute the offset, but \n' +\
                '       there are no initial values available.\n' +\
                '       Please add a model before proceeding.')
          return
       else:
          initial_vals_array = np.array([var['Value'] for var in self.model.initial_variable_list])
          self.nll_extra_offset = 0.
          self.nll_offset = self.ComputeNegLogLikelihood( initial_vals_array )
       return


   ##############################################################################################
   def SetAllVariablesFloating( self ):
       for i in range(len(self.model.variable_list)):
           self.model.variable_list[i]['IsFixed'] = False

   ##############################################################################################
   def SetVariableFixStatus( self, var_name, isFixedInput ):
       var_idx = self.GetVariableIndex( var_name )
       (self.model.variable_list[var_idx])['IsFixed'] = isFixedInput

   ##############################################################################################
   def SetFractionalMinuitInputError( self, var_name, new_minuit_error ):
       # Here, the input should be a fraction (say, 0.05) and we'll
       # scale that to the value of the variable.
       var_idx = self.GetVariableIndex( var_name )
       (self.model.variable_list[var_idx])['MinuitInputError'] =  \
              (self.model.variable_list[var_idx])['Value'] * new_minuit_error

   ##############################################################################################
   def GetVariableFixTuple( self ):
       # iMinuit requires a tuple of booleans to tell it which parameters
       # should be fixed in the fit.
       return tuple( var['IsFixed'] for var in self.model.variable_list )

   ##############################################################################################
   def GetVariableNamesTuple( self ):
       # iMinuit requires a tuple containing the names of each variable
       return tuple( var['Name'] for var in self.model.variable_list )

   ##############################################################################################
   def GetVariableLimitsTuple( self ):
       return tuple( var['Limits'] for var in self.model.variable_list )

   ##############################################################################################
   def GetMinuitInputErrorTuple( self ):
       # iMinuit requires a tuple containing the "error", which
       # is a parameter that I think defines the step size.
       return tuple( var['MinuitInputError'] for var in self.model.variable_list )

   ##############################################################################################
   def SetGaussianConstraintAbsolute( self, var_name, constraint_value, constraint_width ):
       var_idx = self.GetVariableIndex( var_name )
       if not self.model.variable_list[var_idx]['IsConstrained']:
          self.model.variable_list[var_idx]['IsConstrained'] = True
          self.model.constraints.append( {'Name': self.model.variable_list[var_idx]['Name'], \
                                          'Index': var_idx, \
                                          'Value': constraint_value, \
                                          'Width': constraint_width } )
       else:
          for constraint in self.model.constraints:
              if var_name in constraint['Name']:
                 constraint['Value'] = constraint_value
                 constraint['Width'] = constraint_width
                 break
              else:
                 continue                 


   ##############################################################################################
   def SetGaussianConstraintFractional( self, var_name, constraint_value, constraint_width_fractional ):
       var_idx = self.GetVariableIndex( var_name )
       if not self.model.variable_list[var_idx]['IsConstrained']:
          self.model.variable_list[var_idx]['IsConstrained'] = True
          self.model.constraints.append( {'Name': self.model.variable_list[var_idx]['Name'], \
                                          'Index': var_idx, \
                                          'Value': constraint_value, \
                                          'Width': constraint_width_fractional * constraint_value } )
       else:
          for constraint in self.model.constraints:
              if var_name in constraint['Name']:
                 constraint['Value'] = constraint_value
                 constraint['Width'] = constraint_width_fractional *  constraint_value
                 break
              else:
                 continue                 

   ##############################################################################################
   def ClearConstraints( self ):
       self.model.constraints = []
       for var in self.model.variable_list:
           var['IsConstrained'] = False


   ##############################################################################################
   def PlotModelDistributions( self, output_filename='test_plot.png', plot_data=False, save=True, show=False ):

       # Set up the plotting parameters
       plt.rcParams.update({'font.size': 14})
       custom_cycler = cycler( color = [ (1.,0.5,0.),\
                                         (0.,0.,1.,0.5),\
                                         (0.,0.8,0.),\
                                         (1.,0.,0.),\
                                         (0.5,1.,0.5),\
                                         (0.,0.8,0.8),\
                                         (0.1,0.6,0.5),\
                                         (1.,0.,1.),\
                                         (0.,0.,0.,0.2),\
                                         (0.4,0.4,0.4),\
                                         (0.5,0.,0.5) ] )
       plt.rc('axes', prop_cycle=custom_cycler)

       self.fig, self.ax = plt.subplots (2, 2, figsize=(12, 10))

       # Initialize the summed histogram      
       h_sum = hl.hist( [np.array([0.]),np.array([0.]),np.array([0.])] , \
                         bins=self.model.pdfs[0].bins)

       # Loop over pdfs and each to plot
       for i in range(len(self.model.variable_list)):
           var = self.model.variable_list[i]
           if 'Num' in var['Name']:

              weight = var['Value']
              pdf = self.model.pdfs[i]
              component_name = ''.join( var['Name'].split('_')[1:] )
              print('Plotting {}'.format(component_name))

              h_sum += ( weight * pdf )

              hl.plot1d( self.ax[0,0], (weight * pdf)[0:1].project([1]).log10(), label=component_name )
              hl.plot1d( self.ax[0,1], (weight * pdf)[1:].project([1]).log10() )
              hl.plot1d( self.ax[1,1], (weight * pdf)[1:].project([2]).log10() )
              hl.plot1d( self.ax[1,0], (weight * pdf)[0:1].project([2]).log10() )

       hl.plot1d(self.ax[0,0],h_sum[0:1].project([1]).log10(),color='b',label='Total Sum')
       hl.plot1d(self.ax[0,1],h_sum[1:].project([1]).log10(),color='b')
       hl.plot1d(self.ax[1,1],h_sum[1:].project([2]).log10(),color='b')
       hl.plot1d(self.ax[1,0],h_sum[0:1].project([2]).log10(),color='b') 

       if plot_data:
          hl.plot1d( self.ax[0,0], (self.dataset[0:1]).project([1]).log10(),crosses=True,\
                     label='Toy Data Sample',color='k')
          hl.plot1d( self.ax[0,1], (self.dataset[1:]).project([1]).log10(),crosses=True,color='k')
          hl.plot1d( self.ax[1,1], (self.dataset[1:]).project([2]).log10(),crosses=True,color='k')
          hl.plot1d( self.ax[1,0], (self.dataset[0:1]).project([2]).log10(),crosses=True,color='k')

       self.fig.legend(ncol=4,facecolor=(1.,1.,1.),framealpha=1.,loc='upper center')
       self.ax[0,0].set_ylim(-2.,7.)
       self.ax[0,0].set_xlim(700.,3500.)
       self.ax[0,0].set_ylabel('Log10(counts)')
       self.ax[0,0].set_xlabel('Energy (keV)')
       self.ax[0,1].set_ylim(-2.,7.)
       self.ax[0,1].set_xlim(700.,3500.)
       self.ax[0,1].set_ylabel('Log10(counts)')
       self.ax[0,1].set_xlabel('Energy (keV)')
       self.ax[1,1].set_ylim(-2.,7.)
       self.ax[1,1].set_xlim(0.,640.)
       self.ax[1,1].set_ylabel('Log10(counts)')
       self.ax[1,1].set_xlabel('Standoff (mm)')
       self.ax[1,0].set_ylim(-2.,7.)
       self.ax[1,0].set_xlim(0.,640.)
       self.ax[1,0].set_ylabel('Log10(counts)')
       self.ax[1,0].set_xlabel('Standoff (mm)')

       if save:
          plt.savefig(output_filename,dpi=200,bbox_inches='tight')

       if show:
          plt.show()

       return





 
