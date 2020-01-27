import numpy as np
import histlite as hl
import pandas as pd
import nEXOFitModel
import copy

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

   def AddDataset( self, input_dataset ):
       self.dataset = input_dataset
       if self.model is not None:
          self.nll_offset = 0.
          self.SetOffset() 

   def AddPDFDataframeToModel( self, df_pdfs ):
       self.model_obj.AddPDFsFromDataframe( df_pdfs )
       self.model = self.model_obj.GenerateModelDistribution()
       self.variable_list = self.model_obj.variable_list
       for var_row in self.variable_list:
           var_row['IsFixed'] = False
           var_row['MinuitError'] = 0.05*var_row['Value']
       self.initial_values = copy.deepcopy(self.variable_list)

   def ComputeNegLogLikelihood( self, var_values ):
       #
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
                  self.nll_offset 
       return self.nll

   def PrintVariableList( self ):
       print('{:<22} {:<12} {:<9} {:<12}'.format('Variable name:','Value:','IsFixed:','MinuitError:'))
       for var in self.variable_list:
           print('{:<22} {:<12.5} {:<9} {:<12}'.format(var['Name'],var['Value'],str(var['IsFixed']),var['MinuitError']))

   def GetVariableIndex( self, var_name ):
       index = 0
       for i in range(0,len(self.variable_list)):
           if var_name in self.variable_list[i]['Name']:
              index = i
              break
       return index

   def SetOffset( self ):
       if self.initial_values == []:
          print('ERROR: Attempting to compute the offset, but \n' +\
                '       there are no initial values available.\n' +\
                '       Please add a model before proceeding.')
          return
       else:
          initial_vals_array = np.array([var['Value'] for var in self.initial_values])
          self.nll_offset = self.ComputeNegLogLikelihood( initial_vals_array )
       return

   def SetAllVariablesFloating( self ):
       for i in range(len(self.variable_list)):
           self.variable_list[i]['IsFixed'] = False

   def SetVariableFixStatus( self, var_name, isFixedInput ):
       var_idx = self.GetVariableIndex( var_name )
       (self.variable_list[var_idx])['IsFixed'] = isFixedInput

   def SetFractionalMinuitError( self, var_name, new_minuit_error ):
       # Here, the input should be a fraction (say, 0.05) and we'll
       # scale that to the value of the variable.
       var_idx = self.GetVariableIndex( var_name )
       (self.variable_list[var_idx])['MinuitError'] =  \
              (self.variable_list[var_idx])['Value'] * new_minuit_error

   def GetVariableFixTuple( self ):
       # iMinuit requires a tuple of booleans to tell it which parameters
       # should be fixed in the fit.
       return tuple( var['IsFixed'] for var in self.variable_list )

   def GetVariableNamesTuple( self ):
       # iMinuit requires a tuple containing the names of each variable
       return tuple( var['Name'] for var in self.variable_list )

   def GetMinuitErrorTuple( self ):
       # iMinuit requires a tuple containing the "error", which
       # is a parameter that I think defines the step size.
       return tuple( var['MinuitError'] for var in self.variable_list )
