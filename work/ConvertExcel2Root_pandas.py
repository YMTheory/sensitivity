########################################################################
##
## Main script to convert values in summary Excel table into ROOT tree
## Usage: 
## $python ConvertExcel2Root.py
##
########################################################################


import ROOT
import sys
import pandas as pd
import NameDict
import time
import os
import openpyxl

ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')

######################################################
############### EXCEL TABLE READER ###################
######################################################
class ExcelTableReader:
    # This class carries all table info
    # Modify '__init__' appropriately for table in use  

    def __init__(self,inTableName):
        self.filename = inTableName #'../tables/Summary_v68_2016-06-21_0nu.xlsx' #'../tables/Summary_v62_2016-06-04_0nu_tpc_elec.xlsx'
        #writer = pd.ExcelWriter(inTableName,engine='openpyxl')
        #writer.save()
        wb = openpyxl.load_workbook( inTableName )
        wb.template = True
        wb.save( inTableName )

        self.quantile = 1.64

        self.suffixes = ['SS','MS']

        self.specActivSheet = 'SpecificActivities'
        self.countsSheet = '%s_ExpectedCounts'
        self.hitEffSheet = 'MC_RawCounts_%s_Integrals'
        self.halflifeSheet = 'Halflives'                

        print 'Loading sheets...'
        self.dfSpecActiv = self.ReadBasicSheet( self.specActivSheet )
        #self.dfHalflife = self.ReadBasicSheet( self.halflifeSheet )
        self.dfHalflife = pd.read_excel( self.filename, sheet_name=self.halflifeSheet, header=None )
        self.dfHalflife.columns = ['Isotopes','Halflife (yrs)']
        self.dfCountsSS = self.ReadNestedSheet( self.countsSheet % 'SS' )
        self.dfCountsMS = self.ReadNestedSheet( self.countsSheet % 'MS' )
        self.dfHitEffSS = self.ReadNestedSheet( self.hitEffSheet % 'SS' )
        self.dfHitEffMS = self.ReadNestedSheet( self.hitEffSheet % 'MS' )
        #print(self.dfHalflife.head())
        #print(self.dfSpecActiv.head())
        print 'Sheets loaded.'

        self.name_dict = {y:x for x,y in NameDict.NameDict().data.items()}

        self.components = {}
        for index, row in self.dfSpecActiv.iterrows():
            component = row['Component']
            isotope   = row['Isotope']
            mc_id = row['Monte Carlo']
            #pdf       = '%s_%s' % (self.name_dict[component], isotope.replace('-','') )
            pdf       = '%s_%s' % (isotope.replace('-',''),self.name_dict[component])#, isotope )
            self.components[pdf] = ROOT.ExcelTableValues( ROOT.TString(pdf),\
                                                          ROOT.TString(component),\
                                                          ROOT.TString(isotope),\
                                                          ROOT.TString(mc_id) )
            print 'Component: %s\t Isotope: %s\t MC_ID: %s' % (component, isotope, mc_id)


            # Set halflives
            df = self.dfHalflife
            self.components[pdf].SetHalflife( df.loc[ df['Isotopes']==isotope, 'Halflife (yrs)'].iloc[0] )

            # Setting activities.
            df = self.dfSpecActiv
            thisrow = ( df.loc[ df['Component']==component ] ).loc[ df['Isotope']==isotope ]
            specActiv    = thisrow[ 'Specific Activity [mBq/kg]' ].iloc[0]
            specActivErr = thisrow[ 'Error [mBq/kg]' ].iloc[0]
            rawActiv     = thisrow[ 'Activity [Bq]' ].iloc[0]
            rawActivErr  = thisrow[ 'Error [Bq]' ].iloc[0]
            activ_ID     = thisrow[ 'Source' ].iloc[0]
            mc_ID        = thisrow[ 'Monte Carlo' ].iloc[0]
            self.components[pdf].SetActivity( specActiv,\
                                              specActivErr,\
                                              rawActiv,\
                                              rawActivErr,\
                                              ROOT.TString(activ_ID) )                  

            # Setting SS counts
            df = self.dfCountsSS['>700 keV (700, 3500)']['3.648']
            thisrow = (df.loc[ df['Component']==component ]).loc[ df['Isotope']==isotope ]
            self.components[pdf].SetExpectedCounts( thisrow['C.V.'].iloc[0],\
                                                    thisrow['Error'].iloc[0],\
                                                    thisrow['Upper Limit'].iloc[0],\
                                                    'SS' )

            # Setting MS counts
            df = self.dfCountsMS['>700 keV (700, 3500)']['3.648']
            thisrow = (df.loc[ df['Component']==component ]).loc[ df['Isotope']==isotope ]
            self.components[pdf].SetExpectedCounts( thisrow['C.V.'].iloc[0],\
                                                    thisrow['Error'].iloc[0],\
                                                    thisrow['Upper Limit'].iloc[0],\
                                                    'MS' )
           
            # Setting SS hit efficiencies
            # First, set the n and k variables
            df = self.dfHitEffSS['>700 keV (700, 3500)']['3.648']
            thisrow = (df.loc[ df['Component']==component ]).loc[ df['Isotope']==isotope ]
            eff_n = thisrow['No. of Disint.'].iloc[0]
            eff_k = thisrow['C.V.'].iloc[0]
            self.components[pdf].SetHitEfficiency( eff_n, eff_k, 'SS')

            # Next, set the hit efficiencies in each enenrgy ROI
            for energy_range in self.dfHitEffSS.keys():
              if eff_k == 0:
                for i in range(0,4):
                  self.components[pdf].SetHitEfficiencyROI(i,1,1,1,1,1,1,1,1,'SS')
              else:
                ROI_idx = 0
                if 'FWHM' in energy_range:
                   ROI_idx = 0
                if '1-sigma' in energy_range:
                   ROI_idx = 1
                if '2-sigma' in energy_range:
                   ROI_idx = 2
                if '3-sigma' in energy_range:
                   ROI_idx = 3
#                print 'Energy range: %s, ROI: %s' % (energy_range, ROI_idx)
                self.components[pdf].SetHitEfficiencyROI( ROI_idx,\
                                      self.GetHitEfficiencyVal(component, isotope, self.dfHitEffSS[energy_range]['3.648'], eff_k),\
                                      self.GetHitEfficiencyVal(component, isotope, self.dfHitEffSS[energy_range]['3.0'], eff_k),\
                                      self.GetHitEfficiencyVal(component, isotope, self.dfHitEffSS[energy_range]['2.0'], eff_k),\
                                      self.GetHitEfficiencyVal(component, isotope, self.dfHitEffSS[energy_range]['1.0'], eff_k),\
                                      self.GetHitEfficiencyVal(component, isotope, self.dfHitEffSS[energy_range]['3.5'], eff_k),\
                                      self.GetHitEfficiencyVal(component, isotope, self.dfHitEffSS[energy_range]['2.5'], eff_k),\
                                      self.GetHitEfficiencyVal(component, isotope, self.dfHitEffSS[energy_range]['1.5'], eff_k),\
                                      self.GetHitEfficiencyVal(component, isotope, self.dfHitEffSS[energy_range]['0.5'], eff_k),\
                                      'SS' )
                                       
            # Setting MS hit efficiencies
            # First, set the n and k variables
            df = self.dfHitEffMS['>700 keV (700, 3500)']['3.648']
            thisrow = (df.loc[ df['Component']==component ]).loc[ df['Isotope']==isotope ]
            eff_n = thisrow['No. of Disint.'].iloc[0]
            eff_k = thisrow['C.V.'].iloc[0]
            self.components[pdf].SetHitEfficiency( eff_n, eff_k, 'MS')

            # Next, set the hit efficiencies in each enenrgy ROI
            for energy_range in self.dfHitEffMS.keys():
              if eff_k == 0:
                for i in range(0,4):
                  self.components[pdf].SetHitEfficiencyROI(i,1,1,1,1,1,1,1,1,'SS')
              else:
                ROI_idx = 0
                if 'FWHM' in energy_range:
                   ROI_idx = 0
                if '1-sigma' in energy_range:
                   ROI_idx = 1
                if '2-sigma' in energy_range:
                   ROI_idx = 2
                if '3-sigma' in energy_range:
                   ROI_idx = 3
                self.components[pdf].SetHitEfficiencyROI( ROI_idx,\
                                      self.GetHitEfficiencyVal(component, isotope, self.dfHitEffMS[energy_range]['3.648'], eff_k),\
                                      self.GetHitEfficiencyVal(component, isotope, self.dfHitEffMS[energy_range]['3.0'], eff_k),\
                                      self.GetHitEfficiencyVal(component, isotope, self.dfHitEffMS[energy_range]['2.0'], eff_k),\
                                      self.GetHitEfficiencyVal(component, isotope, self.dfHitEffMS[energy_range]['1.0'], eff_k),\
                                      self.GetHitEfficiencyVal(component, isotope, self.dfHitEffMS[energy_range]['3.5'], eff_k),\
                                      self.GetHitEfficiencyVal(component, isotope, self.dfHitEffMS[energy_range]['2.5'], eff_k),\
                                      self.GetHitEfficiencyVal(component, isotope, self.dfHitEffMS[energy_range]['1.5'], eff_k),\
                                      self.GetHitEfficiencyVal(component, isotope, self.dfHitEffMS[energy_range]['0.5'], eff_k),\
                                      'MS' )


            
              
                
    ##########################################################################
    # Reads the sheets from the .xlsx spreadsheet file into pandas dataframes.
    # The sheets with energy/fiducial bins are dicts of dicts of dataframs,
    # allowing us to index them as: data[<energy_bin>][<fid_vol>][C.V.]
    ##########################################################################
    def ReadBasicSheet(self, inSheetName ):
        
        print '\tReading %s...' % inSheetName
        df = pd.read_excel( self.filename, sheet_name = inSheetName  ) 
   
        return df



    def ReadNestedSheet( self, inSheetName ):

        print '\tReading %s' % inSheetName
        header_rows = 4

        # The data in the sheet can be read in directly by skipping the header.

        df = pd.read_excel( self.filename, sheet_name = inSheetName, header = header_rows ) 
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
            if 'Full' in energy_range: continue # This one isn't used
            sheetData[energy_range] = {}
            for fiducial_vol in unique_fiducial_vols:
              df_tmp = pd.DataFrame() # Create an empty dataframe
              for column in global_columns:   # Add in the global columns 
                  df_tmp[column] = df[column]
              for column in local_columns:    # Add in the local columns
                  col_index = self.GetColIndex( energy_range,\
                                                fiducial_vol,\
                                                n_vals_per_fid_vol,\
                                                n_fid_vol_per_en_bin,\
                                                dfheader )
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

        # Using the above, get start index for relevant data values
        jstart = df_subset.loc[ df_subset['Fiducial mass [tonne]'] == fiducial_volume ].index.values[0]

        # Convert the raw index into the pd.read_excel index given to the relevant column
        col_index = int( jstart ) / int( n_vals_per_fid_vol )
        return col_index

    ############################################################################
    # 
    #
    ############################################################################
    def GetHitEfficiencyVal( self, component, isotope, df, eff_k ):
        thisrow = df.loc[ df['Component']==component ].loc[ df['Isotope']==isotope ]
        return thisrow['C.V.'].iloc[0] * 1./eff_k

        

######################################################
############### ROOT TREE WRITER #####################
######################################################
class RootTreeWriter():
    # This class carries info about ROOT tree
    # Adds group information
    # Modify '__init__' appropriately for desired tree

    def __init__(self,outTableName,pdfFilePath):
        self.filename = outTableName #'../tables/Summary_v68_2016-06-21_0nu.root' #'../tables/Summary_v62_2016-06-04_0nu_tpc_elec.root' #'test_new_tree.root'

        self.pdf_filename_pattern ='../histos/Sens_%s.root'
        self.pdf_file_path = pdfFilePath;
    #'/data/data033/exo/software/nEXO_Sensitivity/quick/v5/histos/PNNL/571mm_third/nEXO_Histos_%s.root' # 'individual_histos/nEXO_Histos_%s.root'
        
        self.file = ROOT.TFile.Open(self.filename,'recreate')
        self.tree = ROOT.TTree('ExcelTableValues','Values from Excel Summary Table')
        
        self.table = ROOT.ExcelTableValues()
        
        self.tree.Branch('table',self.table)

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

    def FindGroup(self, pdf):

        for name in self.groups:
            group = self.groups[name]
            for component in group:
                if pdf == component:
                    return name

        print '***** CANNOT FIND GROUP FOR %s *****' % (pdf)
        return None

    def Fill(self, values):
        print 'Filling tree...'

        values.SetGroup(self.FindGroup(values.fPdf))
        #values.SetFileName(self.pdf_filename_pattern % values.fPdf)
        values.SetFileName( self.GeneratePDFFilename(values) )
        #if values.fIsotope == 'bb0n':
        #    values.fFileName = values.fFileName.replace('resol0.005','resol0.010')
        
        self.table.Copy(values)
        self.table.Print()

        self.tree.Fill()

    def Write(self):
        print 'Writing tree into file...'
        self.file.cd()
        self.tree.Write()

    def CloseFile(self):
        print 'Closing file...', self.file.GetName()
        self.file.Close()
  
    def GeneratePDFFilename(self, values):
        mc_id = values.fMC_ID
        files = os.listdir(self.pdf_file_path)
        for f in files:
            if values.fIsotope.replace('-','') in f and mc_id in f:
              #print(self.pdf_file_path + '/' + f)
              return self.pdf_file_path + '/' + f
        

######################################################
############# TABLE TO TREE CONVERTER ################
######################################################
class Excel2RootConverter():
    # This class converts the table into tree

    def __init__(self, excelTable, rootTree):
        self.table = excelTable
        self.tree = rootTree

    def Run(self):
        print 'Converting table to tree...'
        for pdf in self.table.components:
            self.tree.Fill(self.table.components[pdf])
        self.tree.Write()
        self.tree.CloseFile()

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
        print('\tpython ConvertExcel2Root_pandas.py <inputExcelTableName> <outputRootTableName> </path/to/pdf/rootfiles/>')
        sys.exit('\n')
       

    start_time = time.time()    

    excelTable = ExcelTableReader(inTableName) #excelTable.Print()
    rootTree = RootTreeWriter(outTableName,pathToPDFs)

    table2tree = Excel2RootConverter(excelTable,rootTree)
    table2tree.Run()

    end_time = time.time()
    print('Elapsed time = {} seconds ({} minutes).'.format( end_time-start_time, (end_time-start_time)/60. ) )
