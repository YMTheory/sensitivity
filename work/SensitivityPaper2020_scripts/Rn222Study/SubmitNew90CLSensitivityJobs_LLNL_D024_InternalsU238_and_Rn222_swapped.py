#!/usr/local/bin/python

import os
execdir = "/g/g20/lenardo1/nEXO/sensitivity/work/"
#outputdir = "/p/lustre1/lenardo1/sensitivity_output/October28_Rn222_Study_Baseline2019/"
#outputdir = "/p/lustre2/lenardo1/sensitivity_output/Nov20_Rn222_OptimizedBinningV1_RadioassayFluct_D-024/"
outputdir = "/p/lustre2/lenardo1/sensitivity_output/Jan2_Rn222_OptimizedBinningV1_D-024_InternalsU238_and_Rn222_swapped/"
outputname = ""
executable_name = 'Compute90PercentLimit_WilksApprox_Rn222Study_D024_InternalsU238_and_Rn222_swapped.py'
#components_table = '/usr/workspace/wsa/nexo/lenardo1/baseline2019_third_pass/ComponentsTable_D-023_merged-v5_final_cuts.h5' 
#components_table = '/p/vast1/nexo/sensitivity2020/pdfs/component_tables/ComponentsTable_D-023_Optimized_DNN_Standoff_Binning_version1.h5'
components_table = '/p/vast1/nexo/sensitivity2020/pdfs/component_tables/'+\
                'ComponentsTable_D-024_Optimized_DNN_Standoff_Binning_version1_merged-v10b_AllCo60.h5'

rn222_scale_factor = 0.01

iter_num = 1
bkg_shape_err = 0.
num_datasets = 2500
num_jobs = 50
num_datasets_per_job = int(num_datasets/num_jobs)

base = "Run_Rn222_"

for iter_num in range(50,50+num_jobs):
		
		scriptfilename = outputdir +  base + '{:0>4.4}'.format(rn222_scale_factor) + "x_{:03}.sub".format(iter_num)
		os.system( "rm -f " + scriptfilename )
		outfilename = outputdir + base + '{:0>4.4}'.format(rn222_scale_factor) + "x_{:03}.out".format(iter_num)
		os.system( "rm -f " + outfilename )
	
		
		thescript = "#!/bin/bash\n" + \
			"#SBATCH -t 08:00:00\n" + \
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
			"python3 SensitivityPaper2020_scripts/Rn222Study/" + executable_name +\
	                                                " {} {} {} {} {} {} \n".format(\
	                                                iter_num, bkg_shape_err, \
	                                                num_datasets_per_job, components_table, \
                                                        outputdir, rn222_scale_factor) + \
			"export STOPTIME=`date +%s`\n" + \
			"echo Stop time $STOPTIME\n" + \
			"export DT=`expr $STOPTIME - $STARTTIME`\n" + \
			"echo CPU time: $DT seconds\n"
		
		scriptfile = open( scriptfilename, 'w' )
		scriptfile.write( thescript )
		scriptfile.close()
		
		os.system( "sbatch " + scriptfilename )


