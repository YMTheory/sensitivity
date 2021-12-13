#!/usr/local/bin/python

import os

caldate = '21_03_15'
compdate = '21_03_01'
num_its = 200

lam_list = {
	'023': [0],
	'024': [0],
	'025': [0]
}

# lam_list = {
# 	'023': [23, 31, 32, 33, 54, 58, 71, 72, 95, 98, 102, 103],
# 	'024': [12, 13, 16, 18, 20, 23, 24, 28, 30, 35, 36, 39, 40, 43, 45, 46, 47, 48, 49, 50, 52, 53, 54, 55, 56, 58,
# 			59, 60, 62, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 84, 85, 86, 87,
# 			88, 89, 90, 91, 92, 94, 96, 97, 101, 102, 103, 107, 109, 112, 113, 114, 115, 116, 118, 120],
# 	'025': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 16, 17, 18, 19, 22, 23, 24, 25, 27, 28, 30, 32, 34, 37, 38, 39,
# 			40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 52, 55, 56, 57, 58, 60, 62, 63, 64, 65, 66, 67, 70, 71, 72, 73,
# 			75, 82, 83, 84, 85, 86, 87, 88, 89, 92, 93, 95, 99, 101, 102, 103, 105, 108, 109, 116, 117, 119, 121]
# 		}

for database in ['023']: #['023', '024', '025']:
	date = caldate + '_DNN1_' + database

	execdir = "/p/lustre2/czyz1/nexo_sensitivity/work/SensitivityPaper2020_scripts/"
	config_loc = "/p/lustre2/czyz1/nexo_sensitivity/work/config/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1.yaml"
	components_tables = "/p/lustre2/nexouser/czyz1/workdir/components_tables/" + compdate + '/'

	base = "Eres_Crit_Lambda_"
	outputdir = "/p/lustre2/nexouser/czyz1/workdir/lambda/" + date + '/'



	# Number of toy datasets to run for each hypothesis
	num_datasets = 5000
	num_datasets = 25

	for resolution in ['0.008']:
	# for resolution in ['0.008', '0.009', '0.01', '0.011', '0.012', '0.013', '0.014', '0.015', '0.016', '0.017', '0.018']:
		comp_loc = components_tables + "ComponentsTable_D-{}_Energy_Res={}.h5".format(database, resolution)

		for iteration in range(num_its):
			# for num in range(0, 122):
			for num in lam_list[database]:
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
					"python3 /p/lustre2/czyz1/nexo_sensitivity/work/SensitivityPaper2020_scripts/Eres_study/Compute_Crit_Lambda_Eres.py " + \
							"{} ".format(num) + \
							"{} ".format(hyp_val) + \
							"{} ".format(num_datasets) + \
							"{} ".format(outputdir) + \
							"{} ".format(config_loc) + \
							"{} ".format(date) + \
							"{} ".format(comp_loc) + \
							"{} ".format(resolution) + \
							"{} ".format(iteration) + "\n" + \
							"export STOPTIME=`date +%s`\n" + \
					"echo Stop time $STOPTIME\n" + \
					"export DT=`expr $STOPTIME - $STARTTIME`\n" + \
					"echo CPU time: $DT seconds\n"

				scriptfile = open( scriptfilename, 'w' )
				scriptfile.write( thescript )
				scriptfile.close()

				os.system( "sbatch " + scriptfilename )
