import numpy as np
import histlite as hl
import pandas as pd
import nEXOFitModel
import copy
import sys

class nEXOFitLikelihood:

   def __init__( self ):
       # nEXOFitLikelihood is a class that will generate a likelihood function
       # given some arbitrary binned likelihoods and binned PDFs. Both PDFs and
       # datasets are assumed to be in the form of Histlite histograms, and 
       # we assume we're performing a binned maximum likelihood fit.
       self.pdfs = {}
       self.dataset = None
       self.model = None
       self.variable_list = []
       self.initial_values = []
       self.model_obj = nEXOFitModel.nEXOFitModel()
       self.nll = np.nan
       self.nll_offset = 0. # -604307041. # This is a legacy number. Offset should be
                                     # recomputed each time you generate a new dataset. 
       self.nll_extra_offset = 0. # This is used to ensure that the minimum is close to zero
                                  # on the second iteration of the fit.
       self.constraints = []

   ##############################################################################################
   def AddDataset( self, input_dataset ):
       self.dataset = input_dataset
       if self.model is not None:
          self.nll_offset = 0.
          self.SetInitialOffset() 

   ##############################################################################################
   def AddPDFDataframeToModel( self, df_pdfs ):
       self.model_obj.AddPDFsFromDataframe( df_pdfs )
       self.model = self.model_obj.GenerateModelDistribution()
       self.variable_list = self.model_obj.variable_list
       for var_row in self.variable_list:
           var_row['IsFixed'] = False
           var_row['MinuitError'] = 0.05*var_row['Value']
           var_row['IsConstrained'] = False
           var_row['Limits'] = (None, None)
       self.constraints = []
       self.initial_values = copy.deepcopy(self.variable_list)

   ##############################################################################################
   def ComputeNegLogLikelihood( self, var_values ):
     
       fast = False
       self.model_obj.UpdateVariables( var_values )
       self.model = self.model_obj.GenerateModelDistribution(fast=fast)
       if not fast:
          self.model = self.model.values
       self.variable_list = self.model_obj.variable_list
      
       # Here we do a binned negative log-likelihood, assuming
       # each bin is an independent, poisson-distributed variable
       mask = (self.model > 0.)
       self.nll = np.sum( self.model[mask]  - \
                  self.dataset.values[mask] * np.log( self.model[mask] ) ) -\
                  self.nll_offset - self.nll_extra_offset 

       # Here we apply the constraints (if any) to the variables.
       for constraint in self.constraints:
           self.nll += ( var_values[ constraint['Index'] ] - constraint['Value'] )**2 / \
                       ( np.sqrt(2) * constraint['Width']**2 )

       return self.nll

   ##############################################################################################
   def GetVariableValues( self ):
       # Returns a numpy array with the variable values
       array = np.zeros(len(self.variable_list))
       for i in range(0,len(self.variable_list)):
           array[i] = self.variable_list[i]['Value']
       return array

   ##############################################################################################
   def SetVariableLimits( self, var_name, upper_limit=None, lower_limit=None):
       var_idx = self.GetVariableIndex( var_name )
       self.variable_list[var_idx]['Limits'] = (lower_limit,upper_limit)

   ##############################################################################################
   def PrintVariableList( self ):
       print('{:<21} {:<12} {:<9} {:<13} {:<14} {:<13}'.format('Variable name:','Value:',\
                                                        'IsFixed:','MinuitError:','IsConstrained:','Limits:'))
       for var in self.variable_list:
           print('{:<21} {:<12.4} {:<9} {:<13.4} {:<14} ({:4},{:4})'.format(var['Name'], \
                                                                      var['Value'],\
                                                                  str(var['IsFixed']),\
                                                                      var['MinuitError'],\
                                                                  str(var['IsConstrained']),\
                                                                  str(var['Limits'][0]),str(var['Limits'][1]) ))

   ##############################################################################################
   def GetVariableIndex( self, var_name ):
       index = -10000
       for i in range(0,len(self.variable_list)):
           if var_name in self.variable_list[i]['Name']:
              index = i
              break
       if index == -10000:
           #print('\n\nERROR: No variable in the likelihood.variable_list contains the name {}.\n'.format(var_name))
           #print('       Please double-check that the variable names in the ComponentsTable match')
           #print('       the ones you\'re trying to access\n\n')
           err_message =  'No variable in the likelihood.variable_list contains the name {}.\n\n'.format(var_name)
           err_message += '\tPlease double-check that the variable names in the ComponentsTable match\n'
           err_message += '\tthe ones you\'re trying to access\n\n'
           raise ValueError(err_message)

       return index

   ##############################################################################################
   def SetInitialOffset( self ):
       if self.initial_values == []:
          print('ERROR: Attempting to compute the offset, but \n' +\
                '       there are no initial values available.\n' +\
                '       Please add a model before proceeding.')
          return
       else:
          initial_vals_array = np.array([var['Value'] for var in self.initial_values])
          self.nll_extra_offset = 0.
          self.nll_offset = self.ComputeNegLogLikelihood( initial_vals_array )
       return

   ##############################################################################################
   def SetExtraOffset( self, offset ):
       self.nll_extra_offset = offset

   ##############################################################################################
   def SetAllVariablesFloating( self ):
       for i in range(len(self.variable_list)):
           self.variable_list[i]['IsFixed'] = False

   ##############################################################################################
   def SetVariableFixStatus( self, var_name, isFixedInput ):
       var_idx = self.GetVariableIndex( var_name )
       (self.variable_list[var_idx])['IsFixed'] = isFixedInput

   ##############################################################################################
   def SetFractionalMinuitError( self, var_name, new_minuit_error ):
       # Here, the input should be a fraction (say, 0.05) and we'll
       # scale that to the value of the variable.
       var_idx = self.GetVariableIndex( var_name )
       (self.variable_list[var_idx])['MinuitError'] =  \
              (self.variable_list[var_idx])['Value'] * new_minuit_error

   ##############################################################################################
   def GetVariableFixTuple( self ):
       # iMinuit requires a tuple of booleans to tell it which parameters
       # should be fixed in the fit.
       return tuple( var['IsFixed'] for var in self.variable_list )

   ##############################################################################################
   def GetVariableNamesTuple( self ):
       # iMinuit requires a tuple containing the names of each variable
       return tuple( var['Name'] for var in self.variable_list )

   ##############################################################################################
   def GetVariableLimitsTuple( self ):
       return tuple( var['Limits'] for var in self.variable_list )

   ##############################################################################################
   def GetMinuitErrorTuple( self ):
       # iMinuit requires a tuple containing the "error", which
       # is a parameter that I think defines the step size.
       return tuple( var['MinuitError'] for var in self.variable_list )

   ##############################################################################################
   def SetGaussianConstraintAbsolute( self, var_name, constraint_value, constraint_width ):
       var_idx = self.GetVariableIndex( var_name )
       if not self.variable_list[var_idx]['IsConstrained']:
          self.variable_list[var_idx]['IsConstrained'] = True
          self.constraints.append( {'Name': self.variable_list[var_idx]['Name'], \
                                    'Index': var_idx, \
                                    'Value': constraint_value, \
                                    'Width': constraint_width } )
       else:
          for constraint in self.constraints:
              if var_name in constraint['Name']:
                 constraint['Value'] = constraint_value
                 constraint['Width'] = constraint_width
                 break
              else:
                 continue                 


   ##############################################################################################
   def SetGaussianConstraintFractional( self, var_name, constraint_value, constraint_width_fractional ):
       var_idx = self.GetVariableIndex( var_name )
       if not self.variable_list[var_idx]['IsConstrained']:
          self.variable_list[var_idx]['IsConstrained'] = True
          self.constraints.append( {'Name': self.variable_list[var_idx]['Name'], \
                                    'Index': var_idx, \
                                    'Value': constraint_value, \
                                    'Width': constraint_width_fractional * constraint_value } )
       else:
          for constraint in self.constraints:
              if var_name in constraint['Name']:
                 constraint['Value'] = constraint_value
                 constraint['Width'] = constraint_width_fractional *  constraint_value
                 break
              else:
                 continue                 

   ##############################################################################################
   def ClearConstraints( self ):
       self.constraints = []
       for var in self.variable_list:
           var['IsConstrained'] = False
           raise ValueError(message)

       return index

   ##############################################################################################
   def SetInitialOffset( self ):
       if self.initial_values == []:
          print('ERROR: Attempting to compute the offset, but \n' +\
                '       there are no initial values available.\n' +\
                '       Please add a model before proceeding.')
          return
       else:
          initial_vals_array = np.array([var['Value'] for var in self.initial_values])
          self.nll_extra_offset = 0.
          self.nll_offset = self.ComputeNegLogLikelihood( initial_vals_array )
       return

   ##############################################################################################
   def SetExtraOffset( self, offset ):
       self.nll_extra_offset = offset

   ##############################################################################################
   def SetAllVariablesFloating( self ):
       for i in range(len(self.variable_list)):
           self.variable_list[i]['IsFixed'] = False

   ##############################################################################################
   def SetVariableFixStatus( self, var_name, isFixedInput ):
       var_idx = self.GetVariableIndex( var_name )
       (self.variable_list[var_idx])['IsFixed'] = isFixedInput

   ##############################################################################################
   def SetFractionalMinuitError( self, var_name, new_minuit_error ):
       # Here, the input should be a fraction (say, 0.05) and we'll
       # scale that to the value of the variable.
       var_idx = self.GetVariableIndex( var_name )
       (self.variable_list[var_idx])['MinuitError'] =  \
              (self.variable_list[var_idx])['Value'] * new_minuit_error

   ##############################################################################################
   def GetVariableFixTuple( self ):
       # iMinuit requires a tuple of booleans to tell it which parameters
       # should be fixed in the fit.
       return tuple( var['IsFixed'] for var in self.variable_list )

   ##############################################################################################
   def GetVariableNamesTuple( self ):
       # iMinuit requires a tuple containing the names of each variable
       return tuple( var['Name'] for var in self.variable_list )

   ##############################################################################################
   def GetVariableLimitsTuple( self ):
       return tuple( var['Limits'] for var in self.variable_list )

   ##############################################################################################
   def GetMinuitErrorTuple( self ):
       # iMinuit requires a tuple containing the "error", which
       # is a parameter that I think defines the step size.
       return tuple( var['MinuitError'] for var in self.variable_list )

   ##############################################################################################
   def SetGaussianConstraintAbsolute( self, var_name, constraint_value, constraint_width ):
       var_idx = self.GetVariableIndex( var_name )
       if not self.variable_list[var_idx]['IsConstrained']:
          self.variable_list[var_idx]['IsConstrained'] = True
          self.constraints.append( {'Name': self.variable_list[var_idx]['Name'], \
                                    'Index': var_idx, \
                                    'Value': constraint_value, \
                                    'Width': constraint_width } )
       else:
          for constraint in self.constraints:
              if var_name in constraint['Name']:
                 constraint['Value'] = constraint_value
                 constraint['Width'] = constraint_width
                 break
              else:
                 continue                 


   ##############################################################################################
   def SetGaussianConstraintFractional( self, var_name, constraint_value, constraint_width_fractional ):
       var_idx = self.GetVariableIndex( var_name )
       if not self.variable_list[var_idx]['IsConstrained']:
          self.variable_list[var_idx]['IsConstrained'] = True
          self.constraints.append( {'Name': self.variable_list[var_idx]['Name'], \
                                    'Index': var_idx, \
                                    'Value': constraint_value, \
                                    'Width': constraint_width_fractional * constraint_value } )
       else:
          for constraint in self.constraints:
              if var_name in constraint['Name']:
                 constraint['Value'] = constraint_value
                 constraint['Width'] = constraint_width_fractional *  constraint_value
                 break
              else:
                 continue                 

   ##############################################################################################
   def ClearConstraints( self ):
       self.constraints = []
       for var in self.variable_list:
           var['IsConstrained'] = False
