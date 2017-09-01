
import ROOT
libRelPath = '../lib/libnEXOSensitivity.so'
ROOT.gSystem.Load(libRelPath)

import os, time 

jobBody = """
#!/bin/sh
#$ -S /bin/sh
if [ -f /usr/local/sge/boson_nest/common/settings.sh ];then
    source /usr/local/sge/boson_nest/common/settings.sh
fi
source %s
cd %s
python %s %s >& %s
"""

############################
####### RUN OPTIONS ########
############################

all_years = [10.0] #[0.5,1.0,2.5,10.]#[5.0] #[0.5,1.0,2.5,10.]
all_signalCounts = [0.0]

nJobs = 200
nFits = 50  # per job
nRate = 1    # draw mean counts every nRate toy fits
startSeed = 1

setupFile = '/data/data033/exo/software/nEXO_Sensitivity/setup_py27.sh'
progPath = '/data/data033/exo/software/nEXO_Sensitivity/quick/v5/work/'
progName = 'RunSensitivity.py'

progOpts = '-n %s -r %s -y %0.1f -c %0.1f -s %s -d %s -o %s -t %s -m %d --ssfrac-improvement %.1f --rn222-rate-correction %.2f'
baTag = False
scaleBkgd = '' #'10 0'

groupsOff = '' #'Far,FullTpcK40,FullTpcCo60,InternalU238,InternalTh232,LXeBb2n,LXeXe137,LXeRn222,VesselTh232,VesselU238'
ssFracImprovement = 1.0
rn222RateCorrection = 1.0
outDir = '../results'
outName = 'fits_db_v73_2016-09-26_bb2n_0nu_rdm' #'fits_hamamatsu_v61_2016-02-24_Si_Cu_0nu_tpc_cryo_elec_fine_rdm' #  'fits_hamamatsu_v68_2016-06-21_0nu_red8x_fine_rdm' # %(ssFracImprovement) 
treeFileName = '../tables/Summary_v73_2016-09-26_bb2n_0nu.root' # '../tables/Summary_v68_2016-06-21_0nu.root' #'../tables/Summary_v61_2016-02-24_Si_Cu_0nu_tpc_cryo_elec.root' # '../tables/Summary_v68_2016-06-21_0nu_red8x.root' #'../tables/Summary_v62_2016-06-04_0nu_tpc_elec.root' #'../tables/Summary_v61_2016-02-24_Si_Cu_0nu_tpc_resol0.010.root'
jobDir = '../cluster/jobs/'
jobName = '%s_%0.1f_years_%0.1f_counts_%04i'
methodBkgd = ROOT.nEXOSensitivity.kRdmCV #kPosUL

binning = '' #'540 800 3500 126 10 640' #'270 800 3500 31 330 640'

############################
######## MAIN LOOP #########
############################

if baTag:
    progOpts += ' --ba-tag'
if groupsOff != '':
    progOpts += ' --turn-off-groups %s' % (groupsOff)
if scaleBkgd != '':
    progOpts += ' --scale-bkgds %s' % (scaleBkgd)
if  binning != '':
    progOpts += ' --binning %s' % (binning)
    outName += '_bin%s'%(binning.replace(' ','_'))    

for years in all_years:
    for counts in all_signalCounts:
        for jobNo in range(nJobs):
            seed = startSeed + jobNo
            curName = jobName % (outName,years,counts,seed)
            jobFileName = jobDir + curName + '.sh'
            logFileName = jobFileName.replace('.sh','.log')
            
            opts = progOpts % (nFits,nRate,years,counts,seed,outDir,outName,treeFileName,methodBkgd,ssFracImprovement,rn222RateCorrection)
        
            jobFile = open(jobFileName,'w')
            jobFile.write(jobBody % (setupFile,progPath,progName,opts,logFileName))
            jobFile.close()
            
            cmd = 'qsub -q atlas.q -e /dev/null -o /dev/null %s' % (jobFileName)
            print cmd
            os.system(cmd)
            time.sleep(1)
