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
       self.variable_list = []
       self.initial_values = []
       self.model_obj = nEXOFitModel.nEXOFitModel()
       self.nll = np.nan
       self.nll_offset = -604307041.       

   def AddDataset( self, input_dataset ):
       self.dataset = input_dataset

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
