#!/usr/local/bin/python

import os
execdir = "/g/g20/lenardo1/nEXO/sensitivity/work/"
lambda_curve = execdir + "SensitivityPaper2020_notebooks/Critical_lambda_efficiency_parameter_5percent_March19_2020.txt"
outputdir = "/p/lustre1/lenardo1/sensitivity_output/May28_2020_nongaussianity_test/"
outputname = ""
executable_name = 'Compute90PercentLimit_WilksApprox_NonGaussianityStudy.py'
tables_dir = '/g/g20/lenardo1/nEXO/sensitivity/tables/NonGaussianityStudy/' 

iter_num = 1
bkg_shape_err = 1.
num_datasets = 500

base = "Run_"

files_in_tables_dir = os.listdir( tables_dir )
table_files = [thisfile for thisfile in files_in_tables_dir if thisfile.endswith('.h5')]


for thisfile in table_files:

#	print('python3 Compute90PercentLimit_PythonCode.py {} Num_Xe-137 1. {}'.format(num,outputdir))
#	continue

	table_name = thisfile.split('.')[0]
	
	scriptfilename = outputdir +  base + table_name + ".sub"
	os.system( "rm -f " + scriptfilename )
	outfilename = outputdir + base + table_name + ".out"
	os.system( "rm -f " + outfilename )

	#hyp_val = (float(num) + 0.001)/2.
	
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
		"python3 SensitivityPaper2020_scripts/NonGaussianityStudy/" + executable_name +\
                                                " {} {} {} {} {} \n".format(\
                                                iter_num, bkg_shape_err, \
                                                num_datasets, tables_dir + thisfile, outputdir) + \
		"export STOPTIME=`date +%s`\n" + \
		"echo Stop time $STOPTIME\n" + \
		"export DT=`expr $STOPTIME - $STARTTIME`\n" + \
		"echo CPU time: $DT seconds\n"
	
	scriptfile = open( scriptfilename, 'w' )
	scriptfile.write( thescript )
	scriptfile.close()
	
	os.system( "sbatch " + scriptfilename )


