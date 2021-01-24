#!/usr/local/bin/python

import os
execdir = "/g/g20/lenardo1/nEXO/sensitivity/work/"
#outputdir = "/p/lustre2/lenardo1/sensitivity_output/Nov20_Xe137_Study_Baseline2019/"
#outputdir = "/p/lustre2/lenardo1/sensitivity_output//"
outputname = ""
executable_name = 'Compute90PercentLimit_WilksApprox_GammaBackgroundsStudy_D023.py'
#components_table = '/g/g20/lenardo1/nEXO/sensitivity/tables/ComponentsTable_D-005.h5' 
#components_table = '/usr/workspace/nexo/lenardo1/baseline2019_third_pass/ComponentsTable_D-023_merged-v5_final_cuts.h5'
#components_table = '/p/vast1/nexo/sensitivity2020/pdfs/component_tables/ComponentsTable_D-023_Optimized_DNN_Standoff_Binning_version1.h5'
components_table = '/p/vast1/nexo/sensitivity2020/pdfs/component_tables/'+\
			'ComponentsTable_D-023_merged-v10b_Optimized_DNN_Standoff_Binning_version1.h5'

bkg_scale_factor = 0.01

num_datasets = 5000
num_jobs = 100
jobs_offset = 0
num_datasets_per_job = int(num_datasets/num_jobs)

base = "Run_gammaBackgrounds_x"

for iter_num in range(jobs_offset,num_jobs+jobs_offset):
		
		scriptfilename = outputdir +  base + '{:0>4.4}'.format(bkg_scale_factor) + "_{:03}.sub".format(iter_num)
		os.system( "rm -f " + scriptfilename )
		outfilename = outputdir + base + '{:0>4.4}'.format(bkg_scale_factor) + "_{:03}.out".format(iter_num)
		os.system( "rm -f " + outfilename )
	
		
		thescript = "#!/bin/bash\n" + \
			"#SBATCH -t 07:30:00\n" + \
			"#SBATCH -A nuphys\n" + \
			"#SBATCH -e " + outfilename + "\n" + \
			"#SBATCH -o " + outfilename + "\n" + \
			"#SBATCH --mail-type=fail\n" + \
			"#SBATCH -J " + base + "\n" + \
			"#SBATCH --export=ALL \n" + \
			"source /usr/gapps/nexo/setup.sh \n" + \
			"source /g/g20/lenardo1/localpythonpackages/bin/activate \n" + \
			"cd " + execdir + "\n" + \
			"export STARTTIME=`date +%s`\n" + \
			"echo Start time $STARTTIME\n" + \
			"python3 SensitivityPaper2020_scripts/Xe137Study/" + executable_name +\
	                                                " {} {} {} {} {} {} \n".format(\
	                                                iter_num, \
	                                                num_datasets_per_job, components_table, \
                                                        outputdir, bkg_scale_factor) + \
			"export STOPTIME=`date +%s`\n" + \
			"echo Stop time $STOPTIME\n" + \
			"export DT=`expr $STOPTIME - $STARTTIME`\n" + \
			"echo CPU time: $DT seconds\n"
		
		scriptfile = open( scriptfilename, 'w' )
		scriptfile.write( thescript )
		scriptfile.close()
		
		os.system( "sbatch " + scriptfilename )


