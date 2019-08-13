import numpy as np
import histlite as hl
import pandas as pd
import nEXOFitModel

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
       self.model_obj.AddHistogramsFromDataframe( df_pdfs )
       self.model = self.model_obj.GenerateModelDistribution()
       self.variable_list = self.model_obj.variable_list
       self.initial_values = self.variable_list

       return

   def ComputeNegLogLikelihood( self, var_values ):
       #
       self.model_obj.UpdateVariables( var_values )
       self.model = self.model_obj.GenerateModelDistribution()
       self.variable_list = self.model_obj.variable_list
      
       # Here we do a binned negative log-likelihood, assuming
       # each bin is an independent, poisson-distributed variable
       mask = (self.model.values > 0.)
       datmask = (self.dataset.values > 0.)
       self.nll = np.sum( self.model.values[mask]  - \
                  self.dataset.values[mask] * np.log( self.model.values[mask] ) ) -\
                  self.nll_offset 
       return self.nll
