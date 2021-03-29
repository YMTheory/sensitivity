#!/usr/local/bin/python

import os

caldate = '21_03_22'
compdate = '21_03_01'
num_its = 5

for database in ['024']: #['023', '024', '025']:
	date = caldate + '_DNN1_' + database

	execdir = "/p/lustre2/czyz1/nexo_sensitivity/work/SensitivityPaper2020_scripts/"
	config_loc = "/p/lustre2/czyz1/nexo_sensitivity/work/config/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1.yaml"
	components_tables = "/p/lustre2/nexouser/czyz1/workdir/components_tables/" + compdate + '/'

	base = "LiveTime_Crit_Lambda_"
	outputdir = "/p/lustre2/nexouser/czyz1/workdir/lambda/" + date + '/'


	# Number of toy datasets to run for each hypothesis
	num_datasets = 1000
	resolution = '0.008'

	for livetime in ['0.5', '1', '2', '5']:
		comp_loc = components_tables + "ComponentsTable_D-{}_Energy_Res={}.h5".format(database, resolution)

		for iteration in range(num_its):
			for num in range(0, 122):

				basename = base + 'step=' + str(num) + '_livetime=' + livetime + '_num_it=' + str(iteration)

				if not os.path.exists(outputdir):
					os.makedirs(outputdir)
				scriptfilename = outputdir + "Sub/" + basename + ".sub"
				if not os.path.exists(outputdir + "Sub/"):
					os.makedirs(outputdir + "Sub/")
				os.system( "rm -f " + scriptfilename )
				outfilename = outputdir + "Out/" + basename + ".out"
				if not os.path.exists(outputdir + "Out/"):
					os.makedirs(outputdir + "Out/")
				os.system( "rm -f " + outfilename )

				hyp_val = num / 4 + 0.000001

				thescript = "#!/bin/bash\n" + \
					"#SBATCH -t 96:00:00\n" + \
					"#SBATCH -A nuphys\n" + \
					"#SBATCH -e " + outfilename + "\n" + \
					"#SBATCH -o " + outfilename + "\n" + \
					"#SBATCH --mail-type=fail\n" + \
					"#SBATCH -J " + base + "\n" + \
					"#SBATCH -p pbatch \n" + \
					"#SBATCH --export=ALL \n" + \
					"source /usr/gapps/nexo/setup.sh \n" + \
					"source /g/g12/czyz1/nexoenv/bin/activate \n" + \
					"cd " + execdir + "\n" + \
					"export STARTTIME=`date +%s`\n" + \
					"echo Start time $STARTTIME\n" + \
					"python3 /p/lustre2/czyz1/nexo_sensitivity/work/SensitivityPaper2020_scripts/Livetime_study/Compute_Crit_Lambda_LiveTime.py " + \
							"{} ".format(num) + \
							"{} ".format(hyp_val) + \
							"{} ".format(num_datasets) + \
							"{} ".format(outputdir) + \
							"{} ".format(config_loc) + \
							"{} ".format(date) + \
							"{} ".format(comp_loc) + \
							"{} ".format(livetime) + \
							"{} ".format(iteration) + "\n" + \
							"export STOPTIME=`date +%s`\n" + \
					"echo Stop time $STOPTIME\n" + \
					"export DT=`expr $STOPTIME - $STARTTIME`\n" + \
					"echo CPU time: $DT seconds\n"

				scriptfile = open( scriptfilename, 'w' )
				scriptfile.write( thescript )
				scriptfile.close()

				os.system( "sbatch " + scriptfilename )
