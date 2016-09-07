
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

all_years = [0.5,1.0,2.5,5.0,10.]#,12.5,15.0,20.0]
all_signalCounts = [0,5] #2.5,10,20,40,80]#  range(11)

nJobs = 50   # per count, but count = 0 which is 10x more
nFits = 220 # per job
nRate = 1    # draw mean counts every nRate toy fits
startSeed = 1

setupFile = '/data/data033/exo/software/nEXO_Sensitivity/setup_py27.sh'
progPath = '/data/data033/exo/software/nEXO_Sensitivity/quick/v5/work/'
progName = 'RunSensitivity.py'
progOpts = '-n %s -r %s -y %0.1f -c %0.1f -s %s -d %s -o %s -t %s -b '

outDir = '../results'
outName = 'disc_hamamatsu_v62_0nu_eff_tpc_rdm' 
treeFileName = '../tables/Summary_v62_2016-06-04_0nu_tpc.root' #'../tables/Summary_v61_2016-02-24_Si_Cu_0nu_tpc.root'
jobDir = '../cluster/jobs/'
jobName = '%s_%0.1f_years_%0.1f_counts_%04i'

for years in all_years:
    for counts in all_signalCounts:
        totJobs = (nJobs,nJobs*10)[counts == 0]
        for jobNo in range(totJobs):
            seed = startSeed + jobNo
            curName = jobName % (outName,years,counts,seed)
            jobFileName = jobDir + curName + '.sh'
            logFileName = jobFileName.replace('.sh','.log')
            
            opts = progOpts % (nFits,nRate,years,counts,seed,outDir,outName,treeFileName)
        
            jobFile = open(jobFileName,'w')
            jobFile.write(jobBody % (setupFile,progPath,progName,opts,logFileName))
            jobFile.close()
            
            cmd = 'qsub -q atlas.q -e /dev/null -o /dev/null %s' % (jobFileName)
            print cmd
            os.system(cmd)
            time.sleep(1)
