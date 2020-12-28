#!/usr/local/bin/python

import os
execdir = "/g/g20/lenardo1/nEXO/sensitivity/work/SensitivityPaper2020_scripts/BaTagging/"
#outputdir = "/p/lustre1/lenardo1/sensitivity_output/October6_2020_critical_lambda_ba_tagging_no_shape_error_histats_final_cuts/"
outputdir = "/p/lustre2/lenardo1/sensitivity_output/Dec20_2020_CriticalLambda_ba_tagging_background_scaling/"

scaling_2nu = 1.25

base = "BaTagging_Sensitivity_Baseline2019_2nuScaling_{:06.6}_".format(scaling_2nu)

# Number of toy datasets to run for each hypothesis
num_datasets=20000


for num in range(0,400):

	scriptfilename = outputdir + "Sub/" +  base + str(num) + ".sub"
	os.system( "rm -f " + scriptfilename )
	outfilename = outputdir + "Out/" + base + str(num) + ".out"
	os.system( "rm -f " + outfilename )

	if num % 5 == 0:
		hyp_val = (float(num) + 0.0000)/40.
		iteration_num = 0
	elif num % 5 == 1:
		hyp_val = (float(num-1) + 0.0000)/40.
		iteration_num = 1
	elif num % 5 == 2:
		hyp_val = (float(num-2) + 0.0000)/40.
		iteration_num = 2
	elif num % 5 == 3:
		hyp_val = (float(num-3) + 0.0000)/40.
		iteration_num = 3
	elif num % 5 == 4:
		hyp_val = (float(num-4) + 0.0000)/40.
		iteration_num = 4

		
	thescript = "#!/bin/bash\n" + \
		"#SBATCH -t 02:00:00\n" + \
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
		"python3 ComputeCriticalLambdaForFixedNumSignal.py {} {} {} {} {}\n".format( \
				iteration_num, hyp_val, num_datasets, scaling_2nu, outputdir) + \
		"export STOPTIME=`date +%s`\n" + \
		"echo Stop time $STOPTIME\n" + \
		"export DT=`expr $STOPTIME - $STARTTIME`\n" + \
		"echo CPU time: $DT seconds\n"
	
	scriptfile = open( scriptfilename, 'w' )
	scriptfile.write( thescript )
	scriptfile.close()
	
	os.system( "sbatch " + scriptfilename )


