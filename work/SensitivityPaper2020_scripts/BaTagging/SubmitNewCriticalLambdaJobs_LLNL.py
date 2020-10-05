#!/usr/local/bin/python

import os
execdir = "/g/g20/lenardo1/nEXO/sensitivity/work/SensitivityPaper2020_scripts/BaTagging/"
outputdir = "/p/lustre1/lenardo1/sensitivity_output/October4_2020_critical_lambda_ba_tagging/"
outputname = ""

base = "BaTagging_Sensitivity_Baseline2019_"


# Number of toy datasets to run for each hypothesis
num_datasets=5000



for num in range(0,200):

	basename = base + str(num)
	
	scriptfilename = outputdir + "Sub/" +  base + str(num) + ".sub"
	os.system( "rm -f " + scriptfilename )
	outfilename = outputdir + "Out/" + base + str(num) + ".out"
	os.system( "rm -f " + outfilename )

	if num % 2 == 0:
		hyp_val = (float(num) + 0.0001)/8.
		iteration_num = 0
	else:
		hyp_val = (float(num-1) + 0.0001)/8.	
		iteration_num = 1
		
	thescript = "#!/bin/bash\n" + \
		"#SBATCH -t 05:00:00\n" + \
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
		"python3 ComputeCriticalLambdaForFixedNumSignal.py {} {} {} {}\n".format(iteration_num,hyp_val,num_datasets,outputdir) + \
		"export STOPTIME=`date +%s`\n" + \
		"echo Stop time $STOPTIME\n" + \
		"export DT=`expr $STOPTIME - $STARTTIME`\n" + \
		"echo CPU time: $DT seconds\n"
	
	scriptfile = open( scriptfilename, 'w' )
	scriptfile.write( thescript )
	scriptfile.close()
	
	os.system( "sbatch " + scriptfilename )


