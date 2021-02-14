#!/usr/local/bin/python

import os

caldate = '21_01_21'
compdate = '21_01_20'
database = '024'
date = caldate + '_DNN1_' + database

execdir = "/p/lustre2/czyz1/nexo_sensitivity/work/SensitivityPaper2020_scripts/"
config_loc = "/p/lustre2/czyz1/nexo_sensitivity/work/config/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1.yaml"
components_tables = "/p/lustre2/nexouser/czyz1/workdir/components_tables/" + compdate + '/'

base = "Eres_Crit_Lambda_"
outputdir = "/p/lustre2/nexouser/czyz1/workdir/lambda/" + date + '/'


# Number of toy datasets to run for each hypothesis
num_datasets = 5000

# offset = 10

for resolution in ['0.008', '0.011', '0.015', '0.018']:
# for resolution in ['0.008', '0.009', '0.01', '0.011', '0.012', '0.013', '0.014', '0.015', '0.016', '0.017', '0.018']:
	comp_loc = components_tables + "ComponentsTable_D-{}_Energy_Res={}.h5".format(database, resolution)

	for num in range(0, 122):

		basename = base + '_' + str(num)

		if not os.path.exists(outputdir):
			os.makedirs(outputdir)
		scriptfilename = outputdir + date + "/Sub/" + basename + ".sub"
		if not os.path.exists(outputdir + date + "/Sub/"):
			os.makedirs(outputdir + date + "/Sub/")
		os.system( "rm -f " + scriptfilename )
		outfilename = outputdir + date + "/Out/" + basename + ".out"
		if not os.path.exists(outputdir + date + "/Out/"):
			os.makedirs(outputdir + date + "/Out/")
		os.system( "rm -f " + outfilename )

		hyp_val = num / 4

		thescript = "#!/bin/bash\n" + \
			"#SBATCH -t 96:00:00\n" + \
			"#SBATCH -A nuphys\n" + \
			"#SBATCH -e " + outfilename + "\n" + \
			"#SBATCH -o " + outfilename + "\n" + \
			"#SBATCH --mail-type=fail\n" + \
			"#SBATCH -J " + base + "\n" + \
			"#SBATCH -p pdebug \n" + \
			"#SBATCH --export=ALL \n" + \
			"source /usr/gapps/nexo/setup.sh \n" + \
			"source /g/g12/czyz1/nexoenv/bin/activate \n" + \
			"cd " + execdir + "\n" + \
			"export STARTTIME=`date +%s`\n" + \
			"echo Start time $STARTTIME\n" + \
			"python3 /p/lustre2/czyz1/nexo_sensitivity/work/SensitivityPaper2020_scripts/Eres_study/Compute_Crit_Lambda_Eres.py " + \
					"{} ".format(num) + \
					"{} ".format(hyp_val) + \
					"{} ".format(num_datasets) + \
					"{} ".format(outputdir) + \
					"{} ".format(config_loc) + \
					"{} ".format(date) + \
					"{} ".format(comp_loc) + \
					"{} ".format(resolution) + "\n" + \
					"export STOPTIME=`date +%s`\n" + \
			"echo Stop time $STOPTIME\n" + \
			"export DT=`expr $STOPTIME - $STARTTIME`\n" + \
			"echo CPU time: $DT seconds\n"

		scriptfile = open( scriptfilename, 'w' )
		scriptfile.write( thescript )
		scriptfile.close()

		os.system( "sbatch " + scriptfilename )
