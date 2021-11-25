#!/usr/local/bin/python
import sys
import os

sys.path.append('../../../modules')

# res_facts = [(None, 0.008), (0.00412, 0.009), (0.006, .01), (0.00755, .011), (0.00894, .012), (0.01025, .013),
#              (0.01149, .014), (0.01269, .015), (0.01386, .016), (0.015, .017), (0.01612, .018)]
dnn_factors = [0., 0.15, 0.2]
working_dir = "/g/g92/samuele/nEXO/sensitivity/work/SensitivityPaper2020_scripts"
config_file = "../config/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1.yaml"
base_label = 'DNN_factor='
path_to_trees = "/p/vast1/nexo/data/merged-v11-mcid-labels"
date = "21_11_24"
output_dir = "/p/lustre2/nexouser/samuele/histogram_files/" + date

# for (resolution_factor, resolution) in res_facts:
for dnn_factor in dnn_factors:
    label = base_label + str(dnn_factor)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(output_dir + "/Sub"):
        os.makedirs(output_dir + "/Sub")
    if not os.path.exists(output_dir + "/Out"):
        os.makedirs(output_dir + "/Out")
    scriptfilename = output_dir + "/Sub/" + label + ".sub"
    os.system("rm -f " + scriptfilename)

    submit_statement = working_dir + f"/BuildHistogramTableFromTrees.py {config_file} {label} {path_to_trees} {output_dir} "
    if dnn_factor:
        submit_statement += f"-d {dnn_factor}"

    thescript = "#!/bin/bash\n" + \
                "#SBATCH -t 48:00:00\n" + \
                "#SBATCH -A nuphys\n" + \
                "#SBATCH --mail-type=fail\n" + \
                "#SBATCH -J " + label + "\n" + \
                "#SBATCH --export=ALL \n" + \
                "source /usr/gapps/nexo/setup.sh \n" + \
                "source /usr/workspace/samuele/mypyenv/bin/activate \n" + \
                "cd " + working_dir + "\n" + \
                "export STARTTIME=`date +%s`\n" + \
                "echo Start time $STARTTIME\n" + \
                "python3 " + submit_statement + "\n"\
                "export STOPTIME=`date +%s`\n" + \
                "echo Stop time $STOPTIME\n" + \
                "export DT=`expr $STOPTIME - $STARTTIME`\n" + \
                "echo CPU time: $DT seconds\n"

    scriptfile = open(scriptfilename, 'w')
    scriptfile.write(thescript)
    scriptfile.close()

    os.system("sbatch " + scriptfilename)