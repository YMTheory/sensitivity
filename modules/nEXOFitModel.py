import histlite as hl
import pandas as pd
import numpy as np

class nEXOFitModel:

   def __init__( self ):
       
       # Here I'll use lists to ensure that the indices are
       # aligned properly.  
       self.pdfs = []
       self.variable_list = []
       self.full_distribution = None

   def AddPDFsFromDataframe( self, input_df, append=False ):

       self.df_pdfs = input_df
       if not append:
          self.pdfs = []
          self.variable_list = []

       for index, row in input_df.iterrows():
           if (row['Group']=='Total Sum')|\
              (row['Group']=='Off')|\
              (row['Group']=='Far'):
                  continue
           self.pdfs.append( row['Histogram'] )

           this_variable_dict = {}
           this_variable_dict['Name'] = 'Num_{}'.format( row['Group'] ) 
           this_variable_dict['Value'] = row['TotalExpectedCounts']
           self.variable_list.append( this_variable_dict )

   def GenerateModelDistribution( self, fast=False ):
       if len(self.pdfs)==0:
           print('ERROR: no pdfs are loaded into the fitted model.')
           return 

       # Loop over pdfs, weight them, and sum them together.
       distributions = [
           (self.pdfs[i].values if fast else self.pdfs[i])
            * self.variable_list[i]['Value']
           for i in range(len(self.pdfs))
       ]
       self.full_distribution = np.sum(distributions, axis=0)
       return self.full_distribution

   def GenerateDataset( self ):
       if self.full_distribution==None:
          self.GenerateModelDistribution()
 
       # Generates a Poisson sample of each bin
       fake_data_values = np.random.poisson( self.full_distribution.values )   
       # Creates a Histlite histogram with the sampled values and errors     
       self.dataset = hl.Hist( self.full_distribution.bins,\
                                 fake_data_values,\
                                 np.sqrt(fake_data_values) )
       return self.dataset

   def UpdateVariables( self, variable_numpy_array ):

       if not (len(variable_numpy_array)==len(self.variable_list)):
          print('WARNING: Variable array passed to the model object '+\
                'does not match the length of the variable list. '+\
                'Variables will not be updated.\n')     

       for i in range(0,len(variable_numpy_array)):
          self.variable_list[i]['Value'] = np.copy(variable_numpy_array[i])
