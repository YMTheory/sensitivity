#!/usr/local/bin/python

import os
execdir = "/g/g99/jamil1/lustre1/sensitivity/work/SensitivityPaper2020_scripts/"
outputdir = "/g/g99/jamil1/lustre1/sensitivity/work/SensitivityPaper2020_scripts/DiscoveryPotential/test/"
outputname = ""
executable_name = 'Compute90PercentLimit_WilksApprox_DiscoveryPotential.py'
components_table = '/usr/workspace/wsa/nexo/lenardo1/baseline2019_first_pass/ComponentsTable_D-023.h5' 

xe137_scale_factor = 1.0

iter_num = 1
bkg_shape_err = 2.5
num_datasets = 1000
num_jobs = 50
num_datasets_per_job = int(num_datasets/num_jobs)

base = "DiscoveryPotential_NullHypothesis"

for iter_num in range(num_jobs):
		
		scriptfilename = outputdir + 'sub/' + base + "_{:03}.sub".format(iter_num)
		os.system( "rm -f " + scriptfilename )
		outfilename = outputdir + 'out/' + base + "_{:03}.out".format(iter_num)
		os.system( "rm -f " + outfilename )
	
		
		thescript = "#!/bin/bash\n" + \
			"#SBATCH -t 01:00:00\n" + \
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
			"python3 DiscoveryPotential/" + executable_name +\
	                                                " {} {} {} {} {} {} \n".format(\
	                                                iter_num, bkg_shape_err, \
	                                                num_datasets_per_job, components_table, \
                                                        outputdir+'h5/', xe137_scale_factor) + \
			"export STOPTIME=`date +%s`\n" + \
			"echo Stop time $STOPTIME\n" + \
			"export DT=`expr $STOPTIME - $STARTTIME`\n" + \
			"echo CPU time: $DT seconds\n"
		
		scriptfile = open( scriptfilename, 'w' )
		scriptfile.write( thescript )
		scriptfile.close()
		
		os.system( "sbatch " + scriptfilename )


