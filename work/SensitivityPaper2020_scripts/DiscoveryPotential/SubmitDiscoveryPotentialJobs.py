#!/usr/local/bin/python

import os
import numpy as np 
execdir = "/g/g99/jamil1/lustre1/sensitivity/work/SensitivityPaper2020_scripts/"
outputdir = "/g/g99/jamil1/lustre1/sensitivity/work/SensitivityPaper2020_scripts/DiscoveryPotential/Data/version1/"
executable_name = 'Compute_DiscoveryPotential.py'
components_table = '/usr/workspace/wsa/nexo/lenardo1/baseline2019_third_pass/ComponentsTable_D-023_merged-v5_final_cuts.h5' 


subdirs = ['h5','sub','out']
for subdir in subdirs: 
	if not os.path.isdir(outputdir+subdir):
		os.makedirs(outputdir+subdir)


xe137_scale_factor = 1.0

iter_start = 500
bkg_shape_err = 2.5
bb0n_count = 12

livetime = 10

num_datasets = 10000
num_jobs = 500
num_datasets_per_job = int(num_datasets/num_jobs)

base = "DiscoveryPotential_%dCounts_%dyr" % (bb0n_count,livetime)

for iter_num in np.arange(iter_start, iter_start+num_jobs, 1):

	scriptfilename = outputdir + 'sub/' + base + "_{:03}.sub".format(iter_num)
	os.system( "rm -f " + scriptfilename )
	outfilename = outputdir + 'out/' + base + "_{:03}.out".format(iter_num)
	os.system( "rm -f " + outfilename )

	
	thescript = "#!/bin/bash\n" + \
		"#SBATCH -t 00:50:00\n" + \
		"#SBATCH -A nuphys\n" + \
		"#SBATCH -e " + outfilename + "\n" + \
		"#SBATCH -o " + outfilename + "\n" + \
		"#SBATCH --mail-type=fail\n" + \
		"#SBATCH -J " + base + "\n" + \
		"#SBATCH --export=ALL \n" + \
		"source /usr/gapps/nexo/setup.sh \n" + \
		"source /g/g99/jamil1/local/toss_3_x86_64/bin/activate \n" + \
		"cd " + execdir + "\n" + \
		"export STARTTIME=`date +%s`\n" + \
		"echo Start time $STARTTIME\n" + \
		"python3 DiscoveryPotential/" + executable_name + " {} {} {} {} {} {} {} {} \n".format(iter_num, \
																									bkg_shape_err, \
																									num_datasets_per_job, \
																									components_table, \
																									outputdir+'h5/', \
																									xe137_scale_factor, \
																									bb0n_count, \
																									livetime) + \
		"export STOPTIME=`date +%s`\n" + \
		"echo Stop time $STOPTIME\n" + \
		"export DT=`expr $STOPTIME - $STARTTIME`\n" + \
		"echo CPU time: $DT seconds\n"
	
	scriptfile = open( scriptfilename, 'w' )
	scriptfile.write( thescript )
	scriptfile.close()
	
	os.system( "sbatch " + scriptfilename )


