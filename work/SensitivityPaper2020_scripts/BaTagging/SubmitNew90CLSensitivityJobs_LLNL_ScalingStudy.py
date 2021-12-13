#!/usr/local/bin/python

import os
execdir = "/g/g20/lenardo1/nEXO/sensitivity/work/SensitivityPaper2020_scripts/BaTagging/"
lambda_curve = execdir + "CriticalLambdaCurves/critical_lambda_2nu_scaling_study_1e-06.txt"
outputdir = "/p/lustre2/lenardo1/sensitivity_output/Jan1_2021_90CL_ba_tagging_finer_spacing_interp1d_scaling_study/"
outputname = ""

#offset = 0.
print(lambda_curve)
scale_factor = '.'.join((lambda_curve.split('_')[-1]).split('.')[:-1])
#for sub_part in (lambda_curve.split('_')[-1]).split('.')[:-1]:
#   scale_factor = scale_factor + sub_part + '.'
#print(scale_factor)

base = "BaTagging_Sensitivity2020_Mergedv5_90CL_{}_".format(scale_factor)

num_datasets = 500

for num in range(0,50):

	basename = base + str(num)
	
	scriptfilename = outputdir + "Sub/" +  base + str(num) + ".sub"
	os.system( "rm -f " + scriptfilename )
	outfilename = outputdir + "Out/" + base + str(num) + ".out"
	os.system( "rm -f " + outfilename )

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
		"python3 Compute90PercentLimit_PythonCode_ScalingStudy.py {} {} {} {} {} \n".format(\
                                                num, num_datasets, scale_factor, lambda_curve, outputdir) + \
		"export STOPTIME=`date +%s`\n" + \
		"echo Stop time $STOPTIME\n" + \
		"export DT=`expr $STOPTIME - $STARTTIME`\n" + \
		"echo CPU time: $DT seconds\n"
	
	scriptfile = open( scriptfilename, 'w' )
	scriptfile.write( thescript )
	scriptfile.close()
	
	os.system( "sbatch " + scriptfilename )


