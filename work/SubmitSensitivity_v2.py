

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

all_years = [5.0] #[0.5,1.0,2.5,5.0,10.]
all_signalCounts = [0.0]
all_bkgdCounts = [0.25,0.5,1.0,2.0,4.0] #[1.0,2.0,2.5,3.0,4.0]
bkgdGroup = 'InternalTh232' #'VesselTh232' #'VesselTh232' #'Far' #'LXeXe137' #'VesselU238' #'LXeRn222' #'InternalU238'

nJobs = 125
nFits = 100 # per job
nRate = 1 # draw mean counts every nRate toy fits
startSeed = 1

setupFile = '/data/data033/exo/software/nEXO_Sensitivity/setup_py27.sh'
progPath = '/data/data033/exo/software/nEXO_Sensitivity/quick/v2/work/'
progName = 'RunSensitivity.py'
progOpts = '-n %s -m %s -y %0.1f -c %0.1f -r %s -d %s -o %s -t %s -g %s -v %0.2f'

outDir = '../results'
outNamePat = 'fits_hamamatsu_0nu_eff_tpc_rdm_%s_%0.2f' 
treeFileName = '../tables/Summary_v61_2016-02-24_Si_Cu_0nu_tpc.root'
jobDir = '../cluster/jobs/'
jobName = '%s_%0.1f_years_%0.1f_counts_%04i'

for years in all_years:
    for counts in all_signalCounts:
        for bkgd in all_bkgdCounts:
            outName = outNamePat % (bkgdGroup,bkgd)
            for jobNo in range(nJobs):
                seed = startSeed + jobNo
                curName = jobName % (outName,years,counts,seed)
                jobFileName = jobDir + curName + '.sh'
                logFileName = jobFileName.replace('.sh','.log')
                
                opts = progOpts % (nFits,nRate,years,counts,seed,outDir,outName,treeFileName,bkgdGroup,bkgd)
        
                jobFile = open(jobFileName,'w')
                jobFile.write(jobBody % (setupFile,progPath,progName,opts,logFileName))
                jobFile.close()
            
                cmd = 'qsub -q atlas.q -e /dev/null -o /dev/null %s' % (jobFileName)
                print cmd
                os.system(cmd)
                time.sleep(1)
