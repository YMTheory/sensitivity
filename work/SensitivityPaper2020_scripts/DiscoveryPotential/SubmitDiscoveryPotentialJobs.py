#!/usr/local/bin/python

import os
import numpy as np 
import pandas as pd 

version = 17
df = pd.read_csv('dp_versions.csv')
df_version = df[df['Name'].str.contains("Version %d" % version)]

components_table = str(df_version['Components Table'].to_numpy()[0])
yaml_card = str(df_version['YAML'].to_numpy()[0])
outputdir = str(df_version['Data Directory'].to_numpy()[0])

execdir = "/p/lustre1/jamil1/sensitivity/work/SensitivityPaper2020_scripts/"
executable_name = 'Compute_DiscoveryPotential.py'



print(df)
print(df_version)
print(yaml_card)
print(components_table)
print(outputdir)
input()


subdirs = ['h5','sub','out']
for subdir in subdirs: 
	if not os.path.isdir(outputdir+subdir):
		os.makedirs(outputdir+subdir)

iter_start = 0
num_datasets = 10000
num_jobs = 100
num_datasets_per_job = int(num_datasets/num_jobs)

bkg_shape_err = 2.5

# Livetimes = [0.5, 1.0, 2.0, 5.0, 10.0]
# bb0n_counts = {0.5: [0,3,5,8], 
# 			   1.0: [0,4,7,10], 
# 			   2.0: [0,5,8,11], 
# 			   5.0: [0,10,13,16], 
# 			   10.0: [0,15,18,20] 
# 			  } 

Livetimes = [10.0]
bb0n_counts = {10.0: [0,15,18,20] } 
bb0n_counts = {10.0: [0] } 




for livetime in Livetimes[::-1]: 
	for bb0n_count in bb0n_counts[livetime]:
		print(livetime,bb0n_count)
		base = "DiscoveryPotential_%dCounts_%dyr" % (bb0n_count,livetime)

		for iter_num in np.arange(iter_start, iter_start+num_jobs, 1):

			scriptfilename = outputdir + 'sub/' + base + "_{:03}.sub".format(iter_num)
			os.system( "rm -f " + scriptfilename )
			outfilename = outputdir + 'out/' + base + "_{:03}.out".format(iter_num)
			os.system( "rm -f " + outfilename )

			
			thescript = "#!/bin/bash\n" + \
				"#SBATCH -t 3:00:00\n" + \
				"#SBATCH -A nuphys\n" + \
				"#SBATCH -e " + outfilename + "\n" + \
				"#SBATCH -o " + outfilename + "\n" + \
				"#SBATCH --mail-type=fail\n" + \
				"#SBATCH -J " + base + "\n" + \
				"#SBATCH --export=ALL \n" + \
				"source /usr/gapps/nexo/setup.sh \n" + \
				"source /g/g99/jamil1/local/toss_3_x86_64/bin/activate \n" + \
				"cd " + execdir + "\n" + \
				"export STARTTIME=`date +%s`\n" + \
				"echo Start time $STARTTIME\n" + \
				"python3 DiscoveryPotential/" + executable_name + " {} {} {} {} {} {} {} {} \n".format(iter_num, bkg_shape_err, num_datasets_per_job, components_table, outputdir+'h5/', bb0n_count, livetime, yaml_card) + \
				"export STOPTIME=`date +%s`\n" + \
				"echo Stop time $STOPTIME\n" + \
				"export DT=`expr $STOPTIME - $STARTTIME`\n" + \
				"echo CPU time: $DT seconds\n"
			
			scriptfile = open( scriptfilename, 'w' )
			scriptfile.write( thescript )
			scriptfile.close()
			
			os.system( "sbatch " + scriptfilename )


