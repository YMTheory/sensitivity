
import ROOT
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

resolList = [0.005,0.010,0.015,0.025,0.029,0.047] #[0.025,0.030] #[0.017,0.020] #[0.005,0.008,0.010,0.012,0.015] #[r/1000. for r in range(5,16)]
bins = (1,10,20)

setupFile = '/data/data033/exo/software/nEXO_Sensitivity/setup_py27.sh'
progPath = '/data/data033/exo/software/nEXO_Sensitivity/quick/v5/work/'
progName = 'RunSmearHisto.py'
progOpts = '-r %.3f -b %d %d %d -i %s -o %s'


jobDir = '../cluster/jobs/'
jobName = '%s_resol_%.3f'

tree = ROOT.TChain('ExcelTableValues')
tree.Add('../tables/Summary_v61_2016-02-24_Si_Cu_0nu_tpc.root')

for i in range(tree.GetEntries()):
    tree.GetEntry(i)

    pdf = str(tree.fPdf)
    inFile = str(tree.fFileName)
    for resol in resolList:
        outFile = inFile.replace('.root','_resol%.3f.root'%(resol))
        opts = progOpts % ( resol,bins[0],bins[1],bins[2],inFile,outFile)
        
        curName = jobName % (pdf,resol)
        jobFileName = jobDir + curName + '.sh'
        logFileName = jobFileName.replace('.sh','.log')

        jobFile = open(jobFileName,'w')
        jobFile.write(jobBody % (setupFile,progPath,progName,opts,logFileName))
        jobFile.close()

        cmd = 'qsub -q atlas.q -e /dev/null -o /dev/null %s' % (jobFileName)
        print cmd
        os.system(cmd)
        time.sleep(1)
