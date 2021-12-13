#!/usr/local/bin/python

import os
execdir = "/g/g20/lenardo1/nEXO/sensitivity/work/"
outputdir = "/p/lustre2/lenardo1/sensitivity_output/Apr10_SensitivityVsFiducialVolume_D024/"
outputname = "Apr10_FiducialVolumeStudy_merged-v11_OptimizedV1Binning"
executable_name = 'Compute90PercentLimit_WilksApprox_FiducialVolumeStudy_D024.py'
components_table_dir = '/p/vast1/nexo/sensitivity2020/pdfs/component_tables/'
components_table_base = 'ComponentsTable_D-024_merged-v11_'

fiducial_volumes = ['INNER1000kg','INNER1500kg','INNER2000kg','INNER2500kg']



iter_num = 1
bkg_shape_err = 0.
num_datasets = 5000
num_jobs = 100
job_offset = 0
num_datasets_per_job = int(num_datasets/num_jobs)

base = "FidVol_"

for fiducial_volume in fiducial_volumes:
	print('\n\nSubmitting jobs with {} fiducial.\n\n'.format(fiducial_volume))

	components_table = components_table_dir + components_table_base + fiducial_volume + '.h5'

	print('Components table: {}'.format(components_table))

	for iter_num in range(job_offset,job_offset+num_jobs):
			
			scriptfilename = outputdir +  base + '{}'.format(fiducial_volume) + "_{:03}.sub".format(iter_num)
			os.system( "rm -f " + scriptfilename )
			outfilename = outputdir + base + '{}'.format(fiducial_volume) + "x_{:03}.out".format(iter_num)
			os.system( "rm -f " + outfilename )
		
			
			thescript = "#!/bin/bash\n" + \
				"#SBATCH -t 08:00:00\n" + \
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
				"python3 SensitivityPaper2020_scripts/FiducialVolumeStudy/" + executable_name +\
		                                                " {} {} {} {} \n".format(\
		                                                iter_num, \
		                                                num_datasets_per_job, components_table, \
	                                                        outputdir) + \
				"export STOPTIME=`date +%s`\n" + \
				"echo Stop time $STOPTIME\n" + \
				"export DT=`expr $STOPTIME - $STARTTIME`\n" + \
				"echo CPU time: $DT seconds\n"
			
			scriptfile = open( scriptfilename, 'w' )
			scriptfile.write( thescript )
			scriptfile.close()
			
			os.system( "sbatch " + scriptfilename )


