########################################################################
##
## Main script to convert values in summary Excel table into ROOT tree
## Usage: 
## $python ConvertExcel2Root.py
##
########################################################################


import sys
import pandas as pd
import numpy as np
import NameDict
import time
import os
#import openpyxl

######################################################
############### EXCEL TABLE READER ###################
######################################################
class ExcelTableReader:
    # This class carries all table info
    # Modify '__init__' appropriately for table in use  

    def __init__(self,inTableName):

        self.DEBUG=False

        self.filename = inTableName

        self.quantile = 1.64

        #self.suffixes = ['SS','MS']

        self.specActivSheet = 'SpecificActivities'
        self.countsSheet = '%s_ExpectedCounts'
        self.hitEffSheet = 'MC_RawCounts_%s_Integrals'
        self.halflifeSheet = 'Halflives'                

        self.groups = {}

        self.groups['Far'] = [ "U238_OuterCryostatResin",\
                               "U238_OuterCryostatFiber",\
                               "U238_OuterCryostatSupportResin",\
                               "U238_OuterCryostatSupportFiber",\
                               "U238_InnerCryostatResin",\
                               "U238_InnerCryostatFiber",\
                               "U238_InnerCryostatSupportResin",\
                               "U238_InnerCryostatSupportFiber",\
                               "U238_InnerCryostatLiner",\
                               "Th232_OuterCryostatResin",\
                               "Th232_OuterCryostatFiber",\
                               "Th232_OuterCryostatSupportResin",\
                               "Th232_OuterCryostatSupportFiber",\
                               "Th232_InnerCryostatResin",\
                               "Th232_InnerCryostatFiber",\
                               "Th232_InnerCryostatSupportResin",\
                               "Th232_InnerCryostatSupportFiber",\
                               "Th232_InnerCryostatLiner",\
                               "Co60_OuterCryostatResin",\
                               "Co60_OuterCryostatFiber",\
                               "Co60_OuterCryostatSupportResin",\
                               "Co60_OuterCryostatSupportFiber",\
                               "Co60_InnerCryostatResin",\
                               "Co60_InnerCryostatFiber",\
                               "Co60_InnerCryostatSupportResin",\
                               "Co60_InnerCryostatSupportFiber",\
                               "Co60_InnerCryostatLiner",\
                               "K40_OuterCryostatResin",\
                               "K40_OuterCryostatFiber",\
                               "K40_OuterCryostatSupportResin",\
                               "K40_OuterCryostatSupportFiber",\
                               "K40_InnerCryostatResin",\
                               "K40_InnerCryostatFiber",\
                               "K40_InnerCryostatSupportResin",\
                               "K40_InnerCryostatSupportFiber",\
                               "K40_InnerCryostatLiner" ]

        group_vessel = ["HFE",\
                        "TPCVessel",\
                        "TPCSupportCone",\
                        "HVTubes",\
                        "HVCables",\
                        "HVFeedthrough",\
                        "HVFeedthroughCore",\
                        "CalibrationGuideTube1",\
                        "CalibrationGuideTube2" ]

        self.groups['VesselU-238'] = ['U238_%s'%(group_comp) for group_comp in group_vessel]
        self.groups['VesselTh-232'] = ['Th232_%s'%(group_comp) for group_comp in group_vessel]
        
        group_internal = [ "Cathode",\
                           "Bulge",\
                           "FieldRings",\
                           "SupportRodsandSpacers",\
                           "SiPMStaves",\
                           "SiPMElectronics",\
                           "SiPMModuleInterposer",\
                           "SiPMCables",\
                           "SiPMs",\
                           #"ChargeTilesCables",\
                           #"ChargeTilesElectronics",\
                           "ChargeModuleSupport",\
                           "ChargeModuleBacking",\
                           "HVPlunger" ]

        self.groups['InternalU-238'] = ['U238_%s'%(group_comp) for group_comp in group_internal]
        self.groups['InternalTh-232'] = ['Th232_%s'%(group_comp) for group_comp in group_internal]
        
        group_tpc_k40 = ["SupportRodsandSpacers","SiPMModuleInterposer","ChargeTilesBacking"]
        self.groups['FullTpcK-40'] = ['K40_%s'%(group_comp) for group_comp in group_tpc_k40]
        
        group_tpc_co60 = ["ChargeTilesCables","ChargeTilesElectronics","ChargeTilesSupport","ChargeTilesBacking","HVPlunger"]
        self.groups['FullTpcCo-60'] = ['Co60_%s'%(group_comp) for group_comp in group_tpc_co60]
        
        self.groups['ActiveLXeRn-222'] = ["Rn222_ActiveLXe"]
        self.groups['InactiveLXeRn-222'] = ["Rn222_InactiveLXe","Rn222_CathodeRadon"]
        self.groups['InactiveLXeXe-137'] = ["Xe137_InactiveLXe"]
        self.groups['ActiveLXeXe-137'] = ["Xe137_ActiveLXe"]
        self.groups['FullLXeBb2n'] = ["bb2n_FullLXe"]
        self.groups['FullLXeBb0n'] = ["bb0n_FullLXe"]
    
        self.groups['GhostComponents'] = ['K40_%s'%(group_comp) for group_comp in group_internal]
        self.groups['GhostComponents'].extend(['Co60_%s'%(group_comp) for group_comp in group_internal])
        self.groups['GhostComponents'].extend(['K40_%s'%(group_comp) for group_comp in group_vessel])
        self.groups['GhostComponents'].extend(['Co60_%s'%(group_comp) for group_comp in group_vessel])
        self.groups['GhostComponents'].extend(['U238_SolderAnode']) 
        self.groups['GhostComponents'].extend(['Th232_SolderAnode']) 
        self.groups['GhostComponents'].extend(['K40_SolderAnode']) 
        self.groups['GhostComponents'].extend(['Co60_SolderAnode']) 
        self.groups['GhostComponents'].extend(['Ag110m_SolderAnode']) 
        self.groups['GhostComponents'].extend(['U238_SolderSiPM']) 
        self.groups['GhostComponents'].extend(['Th232_SolderSiPM']) 
        self.groups['GhostComponents'].extend(['K40_SolderSiPM']) 
        self.groups['GhostComponents'].extend(['Co60_SolderSiPM']) 
        self.groups['GhostComponents'].extend(['Ag110m_SolderSiPM']) 
        self.groups['GhostComponents'].extend(['B8nu_FullLXe'])
        self.groups['GhostComponents'].extend(['Al26_SupportRodsandSpacers'])
        self.groups['GhostComponents'].extend(['Cs137_FieldRings'])    
        self.groups['GhostComponents'].extend(['U238_ChargeModuleCables'])
        self.groups['GhostComponents'].extend(['Th232_ChargeModuleCables'])
        self.groups['GhostComponents'].extend(['Th232_ChargeModuleElectronics'])
        self.groups['GhostComponents'].extend(['U238_ChargeModuleElectronics'])
        
        # END OF CONSTRUCTOR

    ##########################################################################
    # Reads the necessary data from the .xlsx spreadsheet file into a single
    # pandas dataframe. This replaces the old C++ class.
    ##########################################################################
    def ConvertExcel2DataFrame( self ):

        print('Loading sheets...')
        self.dfSpecActiv = self.ReadBasicSheet( self.specActivSheet )
        #self.dfHalflife = self.ReadBasicSheet( self.halflifeSheet )
        self.dfHalflife = pd.read_excel( self.filename, sheet_name=self.halflifeSheet, header=None )
        self.dfHalflife.columns = ['Isotopes','Halflife (yrs)']
        self.dfCountsSS = self.ReadNestedSheet( self.countsSheet % 'SS' )
        self.dfCountsMS = self.ReadNestedSheet( self.countsSheet % 'MS' )
        self.dfHitEffSS = self.ReadNestedSheet( self.hitEffSheet % 'SS' )
        self.dfHitEffMS = self.ReadNestedSheet( self.hitEffSheet % 'MS' )
        print('Sheets loaded.')

        self.name_dict = {y:x for x,y in NameDict.NameDict().data.items()}

        self.components = pd.DataFrame()
        for index, row in self.dfSpecActiv.iterrows():
            thispdf = pd.Series()

            component = row['Component']
            isotope   = row['Isotope']
            thispdf['PDF'] = '%s_%s' % (isotope.replace('-',''),self.name_dict[component])
            thispdf['Component'] = row['Component']
            thispdf['Isotope'] = row['Isotope']
            thispdf['MC ID'] = row['Monte Carlo']
            print('Component: %s\t Isotope: %s\t MC_ID: %s' % (thispdf['Component'], thispdf['Isotope'], thispdf['MC ID']))


            # Set halflives
            df = self.dfHalflife
            thispdf['Halflife'] = df.loc[ df['Isotopes']==isotope, 'Halflife (yrs)'].iloc[0]

            # Setting activities.
            df = self.dfSpecActiv
            thisrow = ( df.loc[ df['Component']==component ] ).loc[ df['Isotope']==isotope ]
            thispdf['SpecActiv']    = thisrow[ 'Specific Activity [mBq/kg]' ].iloc[0]
            thispdf['SpecActivErr'] = thisrow[ 'Error [mBq/kg]' ].iloc[0]
            thispdf['RawActiv']     = thisrow[ 'Activity [Bq]' ].iloc[0]
            thispdf['RawActivErr']  = thisrow[ 'Error [Bq]' ].iloc[0]
            thispdf['Activity ID']     = thisrow[ 'Source' ].iloc[0]

            # Setting SS counts
            df = self.dfCountsSS['>700 keV (700, 3500)']['3.648']
            thisrowSS = (df.loc[ df['Component']==component ]).loc[ df['Isotope']==isotope ]
            df = self.dfCountsMS['>700 keV (700, 3500)']['3.648']
            thisrowMS = (df.loc[ df['Component']==component ]).loc[ df['Isotope']==isotope ]

            thispdf['Expected Counts'] = thisrowSS['C.V.'].iloc[0] + thisrowMS['C.V.'].iloc[0]
            thispdf['Expected Counts Err'] = np.sqrt( thispdf['Expected Counts'] ) #thisrowSS['Error'].iloc[0] + thisrowMS['Error'].iloc[0]
            thispdf['Expected Counts UL'] = thisrowSS['Upper Limit'].iloc[0] + thisrowMS['Upper Limit'].iloc[0]

           
            # Setting SS hit efficiencies
            # First, set the n and k variables
            df = self.dfHitEffSS['>700 keV (700, 3500)']['3.648']
            thisrowSS = (df.loc[ df['Component']==component ]).loc[ df['Isotope']==isotope ]
            df = self.dfHitEffMS['>700 keV (700, 3500)']['3.648']
            thisrowMS = (df.loc[ df['Component']==component ]).loc[ df['Isotope']==isotope ]
            thispdf['TotalHitEff_N'] = thisrowSS['No. of Disint.'].iloc[0]
            thispdf['TotalHitEff_K'] = thisrowSS['C.V.'].iloc[0] + thisrowMS['C.V.'].iloc[0]

            # Set the group name.
            for group in self.groups:
                if thispdf['PDF'] in self.groups[group]:
                   thispdf['Group'] = group
                   break 

            if self.DEBUG:
              print('Halflife: {}\tSpecActiv: {}\tSpecActivErr: {}\tRawActiv: {}\tAct.ID: {}'.format( \
                    thispdf['Halflife'],\
                    thispdf['SpecActiv'],\
                    thispdf['SpecActivErr'],\
                    thispdf['RawActiv'],\
                    thispdf['Activity ID']))
              print('Expected Cts: {}\tExpected Cts Err: {}\tTotalHitEff: {}'.format( \
                                        thispdf['Expected Counts'],\
                                        thispdf['Expected Counts Err'],\
                                        thispdf['TotalHitEff_K']))



            if self.components.empty:
               self.components = pd.DataFrame(columns=thispdf.index.values)
            self.components.loc[index] = thispdf


            
              
                
    ##########################################################################
    # Reads the sheets from the .xlsx spreadsheet file into pandas dataframes.
    # The sheets with energy/fiducial bins are dicts of dicts of dataframs,
    # allowing us to index them as: data[<energy_bin>][<fid_vol>][C.V.]
    ##########################################################################
    def ReadBasicSheet(self, inSheetName ):
        
        print('\tReading %s...' % inSheetName)
        df = pd.read_excel( self.filename, sheet_name = inSheetName  ) 
   
        return df



    def ReadNestedSheet( self, inSheetName ):

        print('\tReading %s' % inSheetName)
        header_rows = 4

        # The data in the sheet can be read in directly by skipping the header.

        df = pd.read_excel( self.filename, sheet_name = inSheetName, header = header_rows ) 
        #for col in df.columns:
        #  print(col)
        # The header needs some massaging due to the way the 
        # Excel file is currently formatted.
        dfheader_tmp = pd.read_excel( self.filename, sheet_name = inSheetName, skipfooter = df.shape[0]+1 )
        dfheader_tmp = dfheader_tmp.T.reset_index()
        dfheader_tmp.columns = range(0,len(dfheader_tmp.columns))
        dfheader_tmp.columns = dfheader_tmp.loc[ dfheader_tmp[0] == 'Energy Range [keV]' ].iloc[0].tolist()
        dfheader_tmp.columns.name = ''
        dfheader = dfheader_tmp.reset_index(drop=True) # Reset index to start at 0
        index_of_col_headers_row = dfheader.index[ ~dfheader['Fiducial mass [tonne]'].isna() ].tolist()[0]
        index_of_first_data_row = index_of_col_headers_row + 1
        dfheader = dfheader.iloc[index_of_first_data_row:]
        dfheader = dfheader.reset_index(drop=True)  

        energy_ranges = dfheader.loc[ ~dfheader['Energy Range [keV]'].str.contains('Unnamed') ]['Energy Range [keV]'].tolist()
        fiducial_vols = dfheader.loc[ ~dfheader['Fiducial mass [tonne]'].isna() ]['Fiducial mass [tonne]'].tolist()

        unique_energy_ranges = list( set( energy_ranges ) )
        unique_fiducial_vols = list( set( fiducial_vols ) )

        n_vals_per_fid_vol = dfheader.loc[ ~dfheader['Fiducial mass [tonne]'].isna() ].index.values[1] - \
                             dfheader.loc[ ~dfheader['Fiducial mass [tonne]'].isna() ].index.values[0]
        n_fid_vol_per_en_bin = len(unique_fiducial_vols) 

        # Define column names for each dataframe. global_columns are the same in each frame,
        # while the local columns change depending energy and fiducial volume bins.
        global_columns = []
        num_global_columns = 0
        for colname in df.columns: 
            if colname == 'C.V.': break
            global_columns.append(colname) 
            num_global_columns += 1
        local_columns = []
        for colname in df.columns[num_global_columns:num_global_columns+n_vals_per_fid_vol]:
            local_columns.append(colname)

        # Actually fill the data from the sheet into a dict of dicts
        sheetData = {}
        for energy_range in unique_energy_ranges:
            #print('Energy Range: {}'.format(energy_range))
            if 'Full' in energy_range: continue # This one isn't used
            sheetData[energy_range] = {}
            for fiducial_vol in unique_fiducial_vols:
              #print('Fiducial Volume: {}'.format(fiducial_vol))
              df_tmp = pd.DataFrame() # Create an empty dataframe
              for column in global_columns:   # Add in the global columns 
                  #print('Column: {}'.format(column))
                  df_tmp[column] = df[column]
              for column in local_columns:    # Add in the local columns
                  #print('Local Column: {}'.format(column))
                  #print('n_vals_per_fid_vol: {}'.format(n_vals_per_fid_vol))
                  #print('n_fid_vol_per_en_bin: {}'.format(n_fid_vol_per_en_bin))
                  col_index = int( self.GetColIndex( energy_range,\
                                                     fiducial_vol,\
                                                     n_vals_per_fid_vol,\
                                                     n_fid_vol_per_en_bin,\
                                                     dfheader ) )
                  #print('col_index = {}'.format(col_index))
                  if int(col_index) > 0:
                    df_tmp[column] = df['%s.%s' % (column,col_index)]
                  elif int(col_index) == 0:
                    df_tmp[column] = df['%s' % column]
              sheetData[energy_range][str(fiducial_vol)] = df_tmp
                
        return sheetData

    ############################################################################
    # Gets the column name (assigned by pd.read_excel in the header files)
    # for a given fiducial mass and energy range
    ############################################################################
    def GetColIndex( self, energy_range, fiducial_volume, n_vals_per_fid_vol, n_fid_vol_per_energy_bin, dfheader ):

        # Get start and end indices for the given energy bin
        istart = dfheader.loc[ dfheader['Energy Range [keV]'] == energy_range ].index.values[0]
        iend   = istart + n_vals_per_fid_vol * n_fid_vol_per_energy_bin 
        df_subset = dfheader.iloc[istart:iend]

        # Using the above, get start index for relevant data values. If the spreadsheet doesn't contain
        # data for the relevant fiducial cut, return a -1, which will skip this entry in the sheet's dataframe.
        if len( df_subset.loc[ df_subset['Fiducial mass [tonne]'] == fiducial_volume ].index.values ) == 0:
            return -1
        jstart = df_subset.loc[ df_subset['Fiducial mass [tonne]'] == fiducial_volume ].index.values[0]

        # Convert the raw index into the pd.read_excel index given to the relevant column
        col_index = int( jstart ) / int( n_vals_per_fid_vol )
        return col_index
 

######################################################
####################  MAIN  ##########################
######################################################
import sys

if __name__ == "__main__":
   
    
    if len(sys.argv) == 4:
        inTableName = sys.argv[1] # '../tables/Summary_v68_2016-06-21_0nu.xlsx'
        outTableName = sys.argv[2] # '../tables/Summary_v68_2016-06-21_0nu.root'
        pathToPDFs = sys.argv[3]
        if not os.path.exists(pathToPDFs):
           sys.exit('\nERROR: path to PDF files does not exist\n')
    else:
        print('\n\nERROR: ConvertExcel2Root_pandas.py requires 3 arguments');
        print('Usage:')
        print('\tpython ConvertExcel2Root_pandas.py <inputExcelTableName> <outputHDF5FileName> </path/to/pdf/rootfiles/>')
        sys.exit('\n')
       

    start_time = time.time()    

    excelTable = ExcelTableReader(inTableName) #excelTable.Print()
    excelTable.ConvertExcel2DataFrame()

    print( excelTable.components )
    excelTable.components.to_hdf(outTableName,key='Components')

    end_time = time.time()
    print('Elapsed time = {} seconds ({} minutes).'.format( end_time-start_time, (end_time-start_time)/60. ) )
