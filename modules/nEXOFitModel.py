import histlite as hl
import pandas as pd
import numpy as np
import copy

from matplotlib import pyplot as plt


class nEXOFitModel:

   #########################################################################
   def __init__( self ):
       
       # Here I'll use lists to ensure that the indices are
       # aligned properly.  
       self.pdfs = []
       self.variable_list = []
       self.full_distribution = None
       self.constraints = []
       self.signal_efficiency_flag = False
       self.background_shape_var_flag = False
       self.initial_variable_list = []
       self.signal_name = None

   #########################################################################
   def AddPDFsFromDataframe( self, input_df, append=False, replace_existing_variables=True ):

       self.df_pdfs = input_df

       if replace_existing_variables:
           if not append:
              self.pdfs = []
              self.variable_list = []
    
           for index, row in input_df.iterrows():
               if (row['Group']=='Total Sum') | (row['Group']=='Off'):
                      continue
               self.pdfs.append( row['Histogram'] )
    
               this_variable_dict = {}
               this_variable_dict['Name'] = 'Num_{}'.format( row['Group'] ) 
               this_variable_dict['Value'] = row['TotalExpectedCounts']
               self.variable_list.append( this_variable_dict )
           for var in self.variable_list:
               var['IsFixed'] = False
               var['FitError'] = None
               var['MinuitInputError'] = np.sqrt(var['Value'])
               var['IsConstrained'] = False
               var['Limits'] = (None,None)
    
           self.initial_variable_list = copy.deepcopy( self.variable_list )

       else:
           for index, row in input_df.iterrows():
               if (row['Group']=='Total Sum') | (row['Group']=='Off'):
                   continue
               var_idx = self.GetVariableIndexByName( 'Num_{}'.format(row['Group']) )
               self.pdfs[var_idx] = row['Histogram'] 
               self.variable_list[var_idx]['Value'] = row['TotalExpectedCounts'] 

   #########################################################################
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

       # If we're running with the efficiency parameter, use it to scale the signal
       if self.signal_efficiency_flag:
          sig_idx = self.GetVariableIndexByName( 'Bb0n' )
          sig_eff_idx = self.GetVariableIndexByName('Signal_Efficiency')
          distributions[sig_idx] = distributions[sig_idx] * self.variable_list[sig_eff_idx]['Value']

       if self.background_shape_var_flag:
          sig_idx = self.GetVariableIndexByName( 'Bb0n' )
          bkg_shape_idx = self.GetVariableIndexByName('Background_Shape_Error')
          if fast:
             distributions[sig_idx] = distributions[sig_idx] + \
                                      self.variable_list[bkg_shape_idx]['Value'] * self.pdfs[sig_idx].values
          else:
             distributions[sig_idx] = distributions[sig_idx] + \
                                      self.variable_list[bkg_shape_idx]['Value'] * self.pdfs[sig_idx]

       self.full_distribution = np.sum(distributions, axis=0)
       return self.full_distribution


   #########################################################################
   def IncludeSignalEfficiencyVariableInFit( self, include_flag=True ):

       self.signal_efficiency_flag = include_flag

       # GetVariableIndex returns 1000 if the variable is not found 
       # in the variable_list
       try: 
          self.GetVariableIndexByName('Signal_Efficiency')
          # If we need to turn off the signal efficiency parameter, 
          # delete it from the variable list
          if not self.signal_efficiency_flag:
             sig_eff_idx = self.GetVariableIndexByName('Signal_Efficiency')
             del self.variable_list[ sig_eff_idx ]
       except ValueError:
          # If we need to turn on the signal efficiency parameter,
          # add it to the variable list.
          if self.signal_efficiency_flag:
             this_variable_dict = {}
             this_variable_dict['Name'] = 'Signal_Efficiency'
             this_variable_dict['Value'] = 1.
             this_variable_dict['IsFixed'] = False
             this_variable_dict['IsConstrained'] = False
             this_variable_dict['Limits'] = (None,None)
             this_variable_dict['MinuitInputError'] = 0.01
             this_variable_dict['FitError'] = None
             self.variable_list.append(this_variable_dict)
             self.initial_variable_list.append(this_variable_dict)

   #########################################################################
   def IncludeBackgroundShapeVariableInFit( self, include_flag=True ):

       self.background_shape_var_flag = include_flag

       # GetVariableIndex returns 1000 if the variable is not found 
       # in the variable_list
       try: 
          self.GetVariableIndexByName('Background_Shape_Error')
          # If we need to turn off the signal efficiency parameter, 
          # delete it from the variable list
          if not self.background_shape_var_flag:
             sig_eff_idx = self.GetVariableIndexByName('Background_Shape_Error')
             del self.variable_list[ sig_eff_idx ]
       except ValueError:
          # If we need to turn on the signal efficiency parameter,
          # add it to the variable list.
          if self.background_shape_var_flag:
             this_variable_dict = {}
             this_variable_dict['Name'] = 'Background_Shape_Error'
             this_variable_dict['Value'] = 0.
             this_variable_dict['IsFixed'] = False
             this_variable_dict['IsConstrained'] = False
             this_variable_dict['Limits'] = (None,None)
             this_variable_dict['MinuitInputError'] = 0.01
             this_variable_dict['FitError'] = None
             self.variable_list.append(this_variable_dict)
             self.initial_variable_list.append(this_variable_dict)

   #########################################################################
   def GetVariableIndexByName( self, name ):
       index = -10000
       counter = 0
       for i in range(0,len(self.variable_list)):
           if name in self.variable_list[i]['Name']:
              index = i
              counter += 1
       if index == -10000:
           #print('\n\nERROR: No variable in the likelihood.variable_list contains the name {}.\n'.format(var_name))
           #print('       Please double-check that the variable names in the ComponentsTable match')
           #print('       the ones you\'re trying to access\n\n')
           err_message =  'No variable in the likelihood.variable_list contains the name {}.\n\n'.format(name)
           err_message += '\tPlease double-check that the variable names in the ComponentsTable match\n'
           err_message += '\tthe ones you\'re trying to access\n\n'
           raise ValueError(err_message)

       if counter > 1:
          print('WARNING in nEXOFitModel.GetVariableIndexByName():')
          print('\tFound more than one variable matching the string {}'.format(name))

       return index

   #########################################################################
   def GenerateDataset( self ):
       if self.full_distribution is None:
          self.GenerateModelDistribution()
 
       # Generates a Poisson sample of each bin
       fake_data_values = np.random.poisson( self.full_distribution.values )   
       # Creates a Histlite histogram with the sampled values and errors     
       self.dataset = hl.Hist( self.full_distribution.bins,\
                                 fake_data_values,\
                                 np.sqrt(fake_data_values) )
       return self.dataset


   #########################################################################
   def UpdateVariables( self, variable_numpy_array ):

       if not (len(variable_numpy_array)==len(self.variable_list)):
          print('WARNING: Variable array passed to the model object '+\
                'does not match the length of the variable list. '+\
                'Variables will not be updated.\n')     

       for i in range(0,len(variable_numpy_array)):

          self.variable_list[i]['Value'] = np.copy(variable_numpy_array[i])


   
   #########################################################################
   def GetIntegralInBinRange( self, bin_range_numpy_array ):

       num_dimensions = len(bin_range_numpy_array)       

       # bin_range_numpy_array is an array of arrays, with length
       # num_dimensions. Each sub-array contains the indices over
       # which you want to integrate

       if self.full_distribution is None:
          self.GenerateModelDistribution()

       if num_dimensions != len(self.full_distribution.shape):
          print('\nERROR: input array does NOT have the same dimensions as\n')
          print('       this model. Check that your axes are correct.')
          print('       Returning 0.....')
          return 0.
      
       temp_sum = self.full_distribution
 
       for i in range(num_dimensions):
           temp_sum = np.sum( temp_sum[ bin_range_numpy_array[i] ], axis=0)
         
       return temp_sum


   #########################################################################
   def GetComponentIntegralInBinRange( self, var_name, bin_range_numpy_array ):

       num_dimensions = len(bin_range_numpy_array)

       var_idx = self.GetVariableIndexByName( var_name )
       this_pdf = self.pdfs[ var_idx ]
        
       if len(bin_range_numpy_array) != len(this_pdf.values.shape):
          print('\nERROR: input array does NOT have the same dimensions as\n')
          print('       this model. Check that your axes are correct.')
          print('       Returning 0.....')
          return 0.

       temp_sum = this_pdf.values * self.variable_list[var_idx]['Value']
       for i in range(num_dimensions):
           temp_sum = np.sum( temp_sum[ bin_range_numpy_array[i] ], axis=0 )

       return temp_sum





