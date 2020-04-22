#!/usr/local/bin/python

import os
execdir = "/g/g20/lenardo1/nEXO/sensitivity/work/"
lambda_curve = execdir + "SensitivityPaper2020_notebooks/Critical_lambda_efficiency_parameter_5percent_March19_2020.txt"
outputdir = "/p/lustre1/lenardo1/sensitivity_output/March25_2020_90CL_bkg_shape_0p01/"
outputname = ""

bkg_shape_err = 0.01

base = "New_sensitivity_run_test_03_"


for num in range(0,100):

#	print('python3 Compute90PercentLimit_PythonCode.py {} Num_Xe-137 1. {}'.format(num,outputdir))
#	continue
	basename = base + str(num)
	
	scriptfilename = outputdir + "Sub/" +  base + str(num) + ".sub"
	os.system( "rm -f " + scriptfilename )
	outfilename = outputdir + "Out/" + base + str(num) + ".out"
	os.system( "rm -f " + outfilename )

	hyp_val = (float(num) + 0.001)/2.
	
	thescript = "#!/bin/bash\n" + \
		"#SBATCH -t 04:00:00\n" + \
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
		"python3 SensitivityPaper2020_scripts/Compute90PercentLimit_PythonCode.py {} {} {} {} \n".format(\
                                                num, bkg_shape_err, lambda_curve, outputdir) + \
		"export STOPTIME=`date +%s`\n" + \
		"echo Stop time $STOPTIME\n" + \
		"export DT=`expr $STOPTIME - $STARTTIME`\n" + \
		"echo CPU time: $DT seconds\n"
	
	scriptfile = open( scriptfilename, 'w' )
	scriptfile.write( thescript )
	scriptfile.close()
	
	os.system( "sbatch " + scriptfilename )


