#!/usr/local/bin/python

import os
execdir = "/g/g20/lenardo1/nEXO/sensitivity_from_scratch/sensitivity/work/"
outputdir = "/g/g20/lenardo1/nEXO/sensitivity_from_scratch/sensitivity/output/"
bg_table_rootfile = "/g/g20/lenardo1/nEXO/sensitivity_from_scratch/sensitivity/tables/Summary_D-005_v22_2018-11-12.root"
#macro = "/g/g20/lenardo1/Simulations/BACCARAT/TMS/"
outputname = "baseline2017_take2_5000_runs_"
base = "Sensitivity_run_test_02_"


for num in range(0,100):

	basename = base + str(num)
	
	scriptfilename = outputdir + "Sub/" +  base + str(num) + ".sub"
	os.system( "rm -f " + scriptfilename )
	outfilename = outputdir + "Out/" + base + str(num) + ".out"
	os.system( "rm -f " + outfilename )

        hyp_val = (float(num) + 0.001)/2.
	
	thescript = "#!/bin/bash\n" + \
		"#SBATCH -t 24:00:00\n" + \
		"#SBATCH -A ared\n" + \
		"#SBATCH -e " + outfilename + "\n" + \
		"#SBATCH -o " + outfilename + "\n" + \
		"#SBATCH --mail-type=fail\n" + \
		"#SBATCH -J " + base + "\n" + \
		"#SBATCH --export=ALL \n" + \
                "source /usr/gapps/cern/root_v6.08.02/setup \n" + \
                "cd " + execdir + "\n" + \
		"export STARTTIME=`date +%s`\n" + \
		"echo Start time $STARTTIME\n" + \
		"python RunSensitivity.py -n 5000 " + \
                                         "-s 1 " + \
                                         "-b " + \
                                         "-y 10 " + \
                                         "-t " + bg_table_rootfile + " " + \
                                         "-d " + outputdir + " " + \
                                         "-o " + outputname + " " + \
                                         "-c " + str(hyp_val) + " " + \
                                         "--turn-off-groups GhostComponents " + \
                                         "\n" + \
		"export STOPTIME=`date +%s`\n" + \
		"echo Stop time $STOPTIME\n" + \
		"export DT=`expr $STOPTIME - $STARTTIME`\n" + \
		"echo CPU time: $DT seconds\n"
	
	scriptfile = open( scriptfilename, 'w' )
	scriptfile.write( thescript )
	scriptfile.close()
	
	os.system( "sbatch " + scriptfilename )


