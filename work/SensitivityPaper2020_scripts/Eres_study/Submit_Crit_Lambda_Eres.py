#!/usr/local/bin/python

import os
execdir = "/p/lustre2/czyz1/nexo_sensitivty/work/SensitivityPaper2020_scripts/Eres_study/"
outputdir = "/p/lustre2/nexouser/czyz1/output/lambda/"
components_table = "/p/lustre2/nexouser/czyz1/workdir/components_tables"
outputname = ""

base = "Rn222Study_CritLambda_Sensitivity2020_"


# Number of toy datasets to run for each hypothesis
num_datasets=400

rn222_scale_factor = 1.
offset = 10


for num in range(0,1000):

	basename = base + str(num)
	
	scriptfilename = outputdir + "Sub/" +  base + str(num) + ".sub"
	os.system( "rm -f " + scriptfilename )
	outfilename = outputdir + "Out/" + base + str(num) + ".out"
	os.system( "rm -f " + outfilename )

	if num % 10 == 0:
		hyp_val = (float(num) + 0.0000)/40.
		iteration_num = 0 + offset
	elif num % 10 == 1:
		hyp_val = (float(num-1) + 0.0000)/40.
		iteration_num = 1 + offset
	elif num % 10 == 2:
		hyp_val = (float(num-2) + 0.0000)/40.
		iteration_num = 2 + offset
	elif num % 10 == 3:
		hyp_val = (float(num-3) + 0.0000)/40.
		iteration_num = 3 + offset
	elif num % 10 == 4:
		hyp_val = (float(num-4) + 0.0000)/40.
		iteration_num = 4 + offset
	elif num % 10 == 5:
		hyp_val = (float(num-4) + 0.0000)/40.
		iteration_num = 5 + offset
	elif num % 10 == 6:
		hyp_val = (float(num-6) + 0.0000)/40.
		iteration_num = 6 + offset
	elif num % 10 == 7:
		hyp_val = (float(num-7) + 0.0000)/40.
		iteration_num = 7 + offset
	elif num % 10 == 8:
		hyp_val = (float(num-8) + 0.0000)/40.
		iteration_num = 8 + offset
	elif num % 10 == 9:
		hyp_val = (float(num-9) + 0.0000)/40.
		iteration_num = 9 + offset

		
	thescript = "#!/bin/bash\n" + \
		"#SBATCH -t 07:00:00\n" + \
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
		"python3 ComputeCriticalLambdaForFixedNumSignal_Rn222study_D-024.py " + \
				"{} ".format(iteration_num) + \
				"{} ".format(hyp_val) + \
				"{} ".format(num_datasets) + \
				"{} ".format(outputdir) + \
				"{} ".format(components_table) + \
				"{} ".format(rn222_scale_factor) + "\n" + \
		"export STOPTIME=`date +%s`\n" + \
		"echo Stop time $STOPTIME\n" + \
		"export DT=`expr $STOPTIME - $STARTTIME`\n" + \
		"echo CPU time: $DT seconds\n"
	
	scriptfile = open( scriptfilename, 'w' )
	scriptfile.write( thescript )
	scriptfile.close()
	
	os.system( "sbatch " + scriptfilename )


