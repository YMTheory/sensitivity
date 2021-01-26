#!/usr/local/bin/python
import sys
import os
sys.path.append('../../../modules')

res_facts = [( , 0.008), (0.00412, 0.009), (0.006, .01), (0.00755, .011), (0.00894, .012), (0.01025, .013),
			 (0.01149, .014), (0.01269, .015), (0.01386, .016), (0.015, .017), (0.01612, .018)]
working_dir = "/p/lustre2/czyz1/nexo_sensitivity/work/SensitivityPaper2020_scripts"
config_file = "../config/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1.yaml"
base_label = 'Energy_Res='
path_to_trees = "/p/vast1/nexo/data/merged-v10b-mcid-labels"
date = "21_01_20"
output_dir = "/p/lustre2/nexouser/czyz1/workdir/histogram_files/" + date


for (resolution_factor, resolution) in res_facts:
	label = base_label + str(resolution)
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	if not os.path.exists(output_dir + "/Sub"):
		os.makedirs(output_dir + "/Sub")
	if not os.path.exists(output_dir + "/Out"):
		os.makedirs(output_dir + "/Out")
	scriptfilename = output_dir + "/Sub/" + label + ".sub"
	os.system("rm -f " + scriptfilename)


	thescript = "#!/bin/bash\n" + \
		"#SBATCH -t 48:00:00\n" + \
		"#SBATCH -A nuphys\n" + \
		"#SBATCH --mail-type=fail\n" + \
		"#SBATCH -J " + label + "\n" + \
		"#SBATCH --export=ALL \n" + \
		"source /usr/gapps/nexo/setup.sh \n" + \
		"source /g/g12/czyz1/nexoenv/bin/activate \n" + \
		"cd " + working_dir + "\n" + \
		"export STARTTIME=`date +%s`\n" + \
		"echo Start time $STARTTIME\n" + \
		"python3 " + working_dir + "/BuildHistogramTableFromTrees.py {} {} {} {} {}\n".format(\
								config_file, label, path_to_trees, output_dir, resolution_factor) + \
		"export STOPTIME=`date +%s`\n" + \
		"echo Stop time $STOPTIME\n" + \
		"export DT=`expr $STOPTIME - $STARTTIME`\n" + \
		"echo CPU time: $DT seconds\n"
	
	scriptfile = open(scriptfilename, 'w')
	scriptfile.write(thescript)
	scriptfile.close()
	
	os.system("sbatch " + scriptfilename)

