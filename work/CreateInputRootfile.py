########################################################################
##
## Download radioassay values and PDFs directly from the materials 
## database to create a ROOT tree (Modified from ConvertExcel2Root_v9.py)
## 
## Usage: 
## $python CreateInputRootfile.py <detector doc id> <output root file name>
##   <detector doc id>: 
##     D-001: 2015 baseline
##     D-005: 2017 baseline
##
########################################################################
import ROOT
import json, math, time
from db_tools.cabinet import *

ROOT.gSystem.Load('../lib/libnEXOSensitivity.so')

# ROOT
def get_integral(h,en_cut,sd_cut):
  xbinl = h.GetXaxis().FindBin(en_cut[0])
  xbinh = h.GetXaxis().FindBin(en_cut[1])-1
  ybinh = h.GetYaxis().FindBin(sd_cut[0])
  ybinl = h.GetYaxis().FindBin(sd_cut[1])-1
  return h.Integral(xbinl,xbinh,ybinl,ybinh)

def normdist(x,mu,sigma):
  # Cumulative distribution function for the standard normal distribution
  if sigma == 0: return 0
  z = (x-mu)/sigma
  return (1.+math.erf(z/math.sqrt(2.)))/2.

def standardize(rad,iso):
  consts = {
             'U-238': {'ppt': 0.01244, 'ppb': 12.44,  'pg/cm^2': 0.00001244, 'ng/cm^2': 0.01244},
             'Th-232': {'ppt': 0.004073, 'ppb': 4.073,  'pg/cm^2': 0.000004073, 'ng/cm^2': 0.004073},
             'K-40': {'ppt': 0.26532, 'ppb': 265.32,  'pg/cm^2': 0.00026532, 'ng/cm^2': 0.26532},
           }
  for key in consts.keys():
    consts[key]['pg/g'] = consts[key]['ppt']
    consts[key]['ng/g'] = consts[key]['ppb']

  if rad['error_type'] == 'Upper limit (90% C.L.)':
    cv = 0
    err = float(rad['specific_activity'])/1.64
  else:
    cv = float(rad['specific_activity'])
    err = float(rad['error'])

  unit = 'mBq/kg'
  if iso in consts.keys():
    if rad['unit'] in ['ppt', 'pg/g', 'ppb', 'ng/g']:
      cv *= consts[iso][rad['unit']]
      err *= consts[iso][rad['unit']]
      unit = 'mBq/kg'
    elif  rad['unit'] in ['pg/cm^2', 'ng/cm^2']:
      cv *= consts[iso][rad['unit']]
      err *= consts[iso][rad['unit']]
      unit = 'mBq/cm^2'
  return {
          'specific_activity': cv,
          'error': err,
          'error_type': 'Symmetric error (68% C.L.)',
          'altlimit': rad['altlimit'],
          'unit': unit,
          'asym_lower_error': rad['asym_lower_error'],
          'systerr': rad['systerr']
         }


######################################################
############### EXCEL TABLE READER ###################
######################################################
class ExcelTableReader:
    # This class carries all table info
    # Modify '__init__' appropriately for table in use  

    def __init__(self,detector_doc_id):

      self.halflives = {       
          "U-238": 4.468e+09,
          "Th-232": 1.4e+10,
          "K-40": 1.248e+09,
          "Co-60": 5.27,
          "Cs-137": 30.08,
          "Rn-222": 4.468e+09,
          "Bi-214": 4.468e+09,
          "Xe-137": 1e+30,
          "Al-26": 7.17e+05,
          "Ag-110m": 0.6843,
          "bb2n": 2.165e+21,
          "bb0n": 1e+30,
          "B8nu": 1e+30,
      }

      self.en_cuts = {
        'fwhm':  (2428, 2488), 
        '3sig':  (2384, 2532),
        '2sig':  (2408, 2507),
        '1sig':  (2433, 2483),
        'gt700': (700, 3500),
        'full':  (0, 3500),
      }

      self.all_sd_cuts = {
        '2015': {
          '0p5t': (631, 333),
          '1t'  : (631, 256),
          '1p5t': (631, 202),
          '2t'  : (631, 159),
          '2p5t': (631, 122),
          '3t'  : (631,  90),
          '3p5t': (631,  36),
          'fv'  : (631,   0),
        },
        '2017': {
          '0p5t': (558.5, 272.811),
          '1t':   (558.5, 195.860),
          '1p5t': (558.5, 141.864),
          '2t':   (558.5, 98.8712),
          '2p5t': (558.5, 62.5620),
          '3t':   (558.5, 30.8216),
          '3p5t': (558.5, 2.43756),
          'fv':   (558.5, 0),
        },
        '2017v2': {
          '0p5t': (566.65, 278.338),
          '1t':   (566.65, 201.352),
          '1p5t': (566.65, 147.338),
          '2t':   (566.65, 104.334),
          '2p5t': (566.65, 68.0170),
          '3t':   (566.65, 36.2705),
          '3p5t': (566.65, 7.88156),
          'fv':   (566.65, 0),
        },
      }
      geom_tag = '2017v2'
      self.sd_cuts = self.all_sd_cuts[geom_tag]
 
    #  self.ReadContents(detector_doc_id)
    #def ReadContents(self,detector_doc_id):

      detector_doc = get_doc(detector_doc_id) # detector document. 2017 baseline: D-005
      raddocs = {}  # cache for radioassay documents
      mcdocs = {}   # cache for MC documents
 
      # For saving as ROOT format
      self.components = {}

      # Loop over all parts
      for component in detector_doc['components']:
        print component['name'], component['material'],component['montecarloid'], component['radioassayid']

        # Get MC document (One MC doc has the rootfiles for all isotopes associated with this part)
        mainmcdoc = get_doc(component['montecarloid'])

        for isotope in self.halflives.keys():
          pdf = '%s_%s' % (component['name'].replace(' ','').replace('(','').replace(')',''),isotope) # "pdf" is just a string used as key

          # Get raddoc from matDB
          isoname = ''.join(isotope.split('-')).lower()   # Reformat from U-238 to u238
          raddoc = get_measurement(component['radioassayid'],isoname)
          # Add to raddocs only when the spec act exists and non-empty.
          if 'specific_activity' in raddoc.keys() and raddoc['specific_activity'] != '': 
            raddocs[pdf] = raddoc

          # Get PDFs from matDB
          for i,s in enumerate(mainmcdoc['rootfiles']):
            # if the PDF for this component+isotope exists, then download
            if s['isotope'] == isotope:
              #mckey = component['montecarloid']+'_'+isotope      # key used in mcdocs
              dbkey = component['montecarloid']+'.'+str(i+1) # key used in the matdb
              filename = 'nEXO_Histos_%s_%s.root' % (component['name'].replace(' ','').replace('(','').replace(')',''), isotope)  # path to save to
              mcpath = get_mc(dbkey,filename)
              #print 'MC file (%s) for %s is saved to %s' % (pdf, isotope, mcpath)
              mcdoc = s
              if pdf not in mcdocs.keys():
                mcdocs[pdf] = s
              break

          # If both raddoc and mcdoc for this comp+iso are found, then add to self.components
          if pdf in mcdocs.keys() and pdf in raddocs.keys():
            print 'PDF: ', pdf
            self.components[pdf] = ROOT.ExcelTableValues(ROOT.TString(pdf),ROOT.TString(component['name']),ROOT.TString(isotope))
            self.components[pdf].SetHalflife(self.halflives[self.components[pdf].fIsotope])

            # Activities
            stdraddoc = standardize(raddoc,isotope)
            totalmass = float(component['mass'])*float(component['quantity'])
            activity     = float(stdraddoc['specific_activity'])*totalmass
            activity_err = float(stdraddoc['error'])*totalmass
            self.components[pdf].SetActivity(
              float(stdraddoc['specific_activity']), # spec act
              float(stdraddoc['error']),             # spec act err
              activity,                              # act 
              activity_err,                          # act err
              ROOT.TString(component['radioassayid']),         # radioassay matdb id
            )
            
            sens = ROOT.TFile(mcpath,'r')
            for mult in ['SS','MS']:

              # Hit efficiencies
              histo = sens.Get('h_StandoffVsEnergy'+mult+'_Smear')
              hiteff_n = float(mcdoc['numbersimed'])
              hiteff_k = get_integral(histo,self.en_cuts['gt700'],self.sd_cuts['fv'])
              #print 'DEBUG', pdf, hiteff_n, hiteff_k
              ratios = {}
              for en_label in ['fwhm','1sig','2sig','3sig']:
                for sd_label in ['fv','1t','2t','3t','0p5t','1p5t','2p5t','3p5t']:
                  label = en_label+'_'+sd_label
                  ratios[label] = 1 if hiteff_k == 0 else (
                    1.*get_integral(histo,self.en_cuts[en_label],self.sd_cuts[sd_label])/hiteff_k
                  )
                  #print 'DEBUG', pdf, label, ratios[label]
              self.components[pdf].SetHitEfficiency(hiteff_n,hiteff_k,mult)
              self.components[pdf].SetHitEfficiencyROI(0, 
                ratios['fwhm_fv'],ratios['fwhm_3t'],ratios['fwhm_2t'],ratios['fwhm_1t'],
                ratios['fwhm_3p5t'],ratios['fwhm_2p5t'],ratios['fwhm_1p5t'],ratios['fwhm_0p5t'],mult)
              self.components[pdf].SetHitEfficiencyROI(3,
                ratios['3sig_fv'],ratios['3sig_3t'],ratios['3sig_2t'],ratios['3sig_1t'],
                ratios['3sig_3p5t'],ratios['3sig_2p5t'],ratios['3sig_1p5t'],ratios['3sig_0p5t'],mult)
              self.components[pdf].SetHitEfficiencyROI(2,
                ratios['2sig_fv'],ratios['2sig_3t'],ratios['2sig_2t'],ratios['2sig_1t'],
                ratios['2sig_3p5t'],ratios['2sig_2p5t'],ratios['2sig_1p5t'],ratios['2sig_0p5t'],mult)
              self.components[pdf].SetHitEfficiencyROI(1,
                ratios['1sig_fv'],ratios['1sig_3t'],ratios['1sig_2t'],ratios['1sig_1t'],
                ratios['1sig_3p5t'],ratios['1sig_2p5t'],ratios['1sig_1p5t'],ratios['1sig_0p5t'],mult)

              # Expected counts
              hiteff = 1.*hiteff_k/hiteff_n*ratios['fwhm_2t']
              hiteff_err = math.sqrt(hiteff_k*ratios['fwhm_2t']*(hiteff_n-hiteff_k*ratios['fwhm_2t'])/hiteff_n**3)
      
              spy = 31556926.  # seconds in a year
              lifetime = self.components[pdf].fHalflife/math.log(2)
              counts = activity * hiteff * spy
              counts_err = math.sqrt((activity_err*hiteff)**2 + (activity*hiteff_err)**2)*spy
              if lifetime < 1e6: 
                corr = lifetime*(1.-math.exp(-1./lifetime))
                counts *= corr
                counts_err *= corr

              G0 = normdist(0,counts,counts_err)
              l_up_2 = ROOT.Math.normal_quantile(0.5*(1+0.9*(1-G0)),counts_err)
              l_up_1 = ROOT.Math.normal_quantile(0.9+(1-0.9)*G0,counts_err)
              l_low = 2*counts-l_up_2
              lowerlimit = 0 if counts_err <= 0 else max(0,l_low)
              upperlimit = 0 if counts_err <= 0 else (1.64*counts_err if counts < 0 else max(l_up_1,l_up_2))
              self.components[pdf].SetExpectedCounts(counts, counts_err, upperlimit, mult)
              # For details on the calculation of lower and upper limits, 
              # see slide 3 of https://ntpc.ucllnl.org/nEXO/images/9/91/160204-errors-radio.pdf 

      return        

    def Print(self):

        for pdf in self.components:
            self.components[pdf].Print()

        return

######################################################
############### ROOT TREE WRITER #####################
######################################################
class RootTreeWriter():
    # This class carries info about ROOT tree
    # Adds group information
    # Modify '__init__' appropriately for desired tree

    def __init__(self,outTableName):
        self.filename = outTableName

        #self.pdf_filename_pattern = '/home/tsan630/Code/user_projects/rhmtsang/MaterialDatabase/CabinetAPI/nEXO_Histos_%s.root'
        self.pdf_filename_pattern = './nEXO_Histos_%s.root'
        
        self.file = ROOT.TFile.Open(self.filename,'recreate')
        self.tree = ROOT.TTree('ExcelTableValues','Values from Excel Summary Table')
        
        self.table = ROOT.ExcelTableValues()
        
        self.tree.Branch('table',self.table)

    def FindGroup(self, pdf):

        # If this pdf involves a special volume (e.g. LXe) or a special isotope (e.g. Rn-222)
        special_pdfs = {
          'ActiveLXe_Xe-137': 'ActiveLXeXe-137',
          'ActiveLXe_Rn-222': 'ActiveLXeRn-222',
          'InactiveLXe_Xe-137': 'InactiveLXeXe-137',
          'InactiveLXe_Rn-222':    'InactiveLXeRn-222',
          'Cathode(Radon)_Rn-222': 'InactiveLXeRn-222',
          'FullLXe_bb2n': 'LXeBb2n',
          'FullLXe_bb0n': 'LXeBb0n',
          'FullLXe_B8nu': 'LXeB8nu',
          'Solder(Anode)_Ag-110m': 'InternalAg-110m',
          'Solder(SiPM)_Ag-110m':  'InternalAg-110m',
          'FieldRings_Cs-137': 'GhostComponents',
          'ActiveLXe_Kr-85': 'GhostComponents',
          'InactiveLXe_Kr-85': 'GhostComponents',
        }

        if pdf in special_pdfs.keys():
          return special_pdfs[pdf]

        # If this is just a regular isotope (U/Th/K/Co)
        component, isotope = pdf.split('_')

        group_far = ['OuterCryostat(Fiber)', 'OuterCryostat(Resin)', 
                     'OuterCryostatSupport(Fiber)', 'OuterCryostatSupport(Resin)',
                     'InnerCryostat(Fiber)', 'InnerCryostat(Resin)', 
                     'InnerCryostatSupport(Fiber)', 'InnerCryostatSupport(Resin)',
                     'InnerCryostatLiner', ]

        group_internal = ['Cathode', 'Bulge', 'SupportRodsandSpacers', 'FieldRings', 
                          'ChargeTilesCables', 'ChargeTilesElectronics', 'ChargeTilesBacking', 'ChargeTilesSupport',
                          'SiPMCables', 'SiPMElectronics', 'SiPMs', 'SiPMModule(Interposer)', 'SiPMStaves',
                          'HVPlunger', 'Solder(Anode)', 'Solder(SiPM)',]

        group_vessel = [ 'HVCables', 'HVFeedthroughCore', 'HVFeedthrough', 'HVTubes',
                         'CalibrationGuideTube1', 'CalibrationGuideTube2', 'TPCSupportCone', 'TPCVessel', 'HFE',]

        if component in group_far: return 'Far'

        if component in group_internal: 
          if isotope == 'U-238':  return 'InternalU-238'
          if isotope == 'Th-232': return 'InternalTh-232'
          if isotope == 'K-40':   return 'FullTpcK-40'
          if isotope == 'Co-60':  return 'FullTpcCo-60'

        if component in group_vessel: 
          if isotope == 'U-238':  return 'VesselU-238'
          if isotope == 'Th-232': return 'VesselTh-232'
          if isotope == 'K-40':   return 'FullTpcK-40'
          if isotope == 'Co-60':  return 'FullTpcCo-60'

        print '***** CANNOT FIND GROUP FOR %s *****' % (pdf)
        return None

    def Fill(self, values):
        print 'Filling tree...'

        values.SetGroup(self.FindGroup(values.fPdf))
        values.SetFileName(self.pdf_filename_pattern % values.fPdf)
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
   
  if len(sys.argv) == 3:
    detector_doc_id = sys.argv[1] # D-005
    outTableName = sys.argv[2] # '../tables/Summary_v68_2016-06-21_0nu.root'
  else:
    print 'Usage: python %s <detector_doc_id> <output root file name>' % sys.argv[0]
    sys.exit(0)

  excelTable = ExcelTableReader(detector_doc_id)
  rootTree = RootTreeWriter(outTableName)
  
  table2tree = Excel2RootConverter(excelTable,rootTree)
  table2tree.Run()
