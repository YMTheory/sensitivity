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
       print('{:<25} {}'.format('Variable name:','Value:'))
       for var in self.variable_list:
           print('{:<25} {:4.4}'.format(var['Name'],var['Value']))

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
