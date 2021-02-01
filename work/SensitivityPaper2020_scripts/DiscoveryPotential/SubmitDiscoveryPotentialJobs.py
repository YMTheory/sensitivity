#!/usr/local/bin/python

import os
import numpy as np 
import pandas as pd 

CLUSTER='SLAC'
version = 39

df = pd.read_csv('dp_versions.csv')
df_version = df[df['Name'].str.contains("Version %d" % version)]

components_table = str(df_version['Components Table'].to_numpy()[0])
yaml_card = str(df_version['YAML'].to_numpy()[0])
outputdir = str(df_version['Data Directory'].to_numpy()[0])

execdir = "/p/lustre1/jamil1/sensitivity/work/SensitivityPaper2020_scripts/"
executable_name = 'Compute_DiscoveryPotential.py'

print(df_version)
print(yaml_card)
print(components_table)
print(outputdir)

# review inputs before submitting jobs
input()

subdirs = ['h5','sub','out']
for subdir in subdirs: 
	if not os.path.isdir(outputdir+subdir):
		os.makedirs(outputdir+subdir)

bkg_shape_err = 2.5
iter_start = 0

# The null hypothesis will contain 100k toy datasets distributed over 1000 jobs with 100 toy data sets per job 
num_datasets = {0:100000}
num_jobs = {0:1000}
num_datasets_per_job = {0:100}

# Each alternative hypothesis will contain just 3k toy datasets distributed over 100 jobs with 30 toy data sets per job 
for ii in range(1,100,1):
	num_datasets[ii] = 5000
	num_jobs[ii] = 50
	num_datasets_per_job[ii] = int(num_datasets[ii]/num_jobs[ii])


Livetimes = [0.5, 1.0, 2.0, 5.0, 10.0]
bb0n_counts = {0.5: [0,1,2,3,4,5,6,7,8],
			   1.0: [0,2,4,6,8,10],
			   2.0: [0,4,6,8,10,12,14],
			   5.0: [0,6,8,10,12,14,16],
			   10.0: [0,10,12,14,16,18,20]} 
# bb0n_counts = {10.0: [0]} 

Livetimes = [10.0]


for livetime in Livetimes[::-1]: 
	for bb0n_count in bb0n_counts[livetime][::-1]:
		print(livetime, bb0n_count, num_jobs[bb0n_count], num_datasets_per_job[bb0n_count])

		base = "DiscoveryPotential_%dCounts_%dyr" % (bb0n_count,livetime)

		for iter_num in np.arange(iter_start, iter_start + num_jobs[bb0n_count], 1):

			scriptfilename = outputdir + 'sub/' + base + "_{:03}.sub".format(iter_num)
			outfilename = outputdir + 'out/' + base + "_{:03}.out".format(iter_num)
			os.system( "rm -f " + scriptfilename )
			os.system( "rm -f " + outfilename )

			if CLUSTER=='LLNL':
				thescript = "#!/bin/bash\n" + \
					"#SBATCH -t 3:00:00\n" + \
					"#SBATCH -A nuphys\n" + \
					"#SBATCH -e " + outfilename + "\n" + \
					"#SBATCH -o " + outfilename + "\n" + \
					"#SBATCH --mail-type=fail\n" + \
					"#SBATCH -J " + base + "\n" + \
					"#SBATCH --export=ALL \n" + \
					"source /usr/gapps/nexo/setup.sh \n" + \
					"conda activate sensitivity2020 \n" + \
					"cd " + execdir + "\n" + \
					"export STARTTIME=`date +%s`\n" + \
					"echo Start time $STARTTIME\n" + \
					"python3 DiscoveryPotential/" + executable_name + " {} {} {} {} {} {} {} {} \n".format(iter_num, bkg_shape_err, num_datasets_per_job[bb0n_count], components_table, outputdir+'h5/', bb0n_count, livetime, yaml_card) + \
					"export STOPTIME=`date +%s`\n" + \
					"echo Stop time $STOPTIME\n" + \
					"export DT=`expr $STOPTIME - $STARTTIME`\n" + \
					"echo CPU time: $DT seconds\n"
	
				scriptfile = open(scriptfilename, 'w')
				scriptfile.write(thescript)
				scriptfile.close()
				os.system("sbatch " + scriptfilename)
			
			elif CLUSTER=='SLAC':
				os.system( "bsub -W 100 -n 1 -R rhel60 -o {} python Compute_DiscoveryPotential.py ".format(outfilename) + \
							" {} {} {} {} {} {} {} {}".format(iter_num, bkg_shape_err, num_datasets_per_job[bb0n_count], components_table, outputdir+'h5/', bb0n_count, livetime, yaml_card))


