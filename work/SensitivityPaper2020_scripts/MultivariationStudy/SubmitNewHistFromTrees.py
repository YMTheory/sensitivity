#!/usr/local/bin/python
import sys
import os

sys.path.append('../../../modules')

# res_factors = [(None, 0.008), (0.00412, 0.009), (0.006, .01), (0.00755, .011), (0.00894, .012), (0.01025, .013),
#               (0.01149, .014), (0.01269, .015), (0.01386, .016), (0.015, .017), (0.01612, .018)]
# dnn_factors = [0., 0.15, 0.177, 0.2, 0.25]

res_factors = [(None, 0.008), (0.006, .01), (0.00755, .011), (0.00894, .012), (0.01149, .014)]
dnn_factors = [0., 0.15, 0.177, 0.2, 0.25]
working_dir = "/g/g92/samuele/nEXO/sensitivity/work/SensitivityPaper2020_scripts"

#config_file = working_dir + "MultivariationStudy/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1.yaml"
# fine_binning = False

# Uncomment the following lines for finer binning, 
# which is useful for getting correct ROI backgrounds (not for sensitivity fits)
config_file = working_dir + "/MultivariationStudy/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1_fineBinning.yaml"
fine_binning = True

path_to_trees = "/p/vast1/nexo/data/merged-v11-mcid-labels"
date = "22_03_30"
output_dir = "/p/lustre2/nexouser/samuele/histogram_files/" + date

# Note: the DNN smoothing takes about ~2 hours just by itself for all files
for (resolution_factor, resolution) in res_factors:
    for dnn_factor in dnn_factors:
        label = f'DNN_factor={dnn_factor}_ERes={resolution}'
        if fine_binning:
            label = label + '_fineBinning'

        os.makedirs(f'{output_dir}/logs/', exist_ok=True)
        os.makedirs(f'{output_dir}/jobs/', exist_ok=True)
        os.makedirs(f'{output_dir}/output/', exist_ok=True)

        scriptfilename = output_dir + "/jobs/" + label + ".sh"
        os.system("rm -f " + scriptfilename)

        submit_statement = working_dir + f"/BuildHistogramTableFromTrees.py {config_file} {label} {path_to_trees} {output_dir}/output/ "
        if dnn_factor:
            submit_statement += f"-d {dnn_factor} "
        if resolution_factor:
            submit_statement += f"-r {resolution_factor} "

        thescript = "#!/bin/bash\n" + \
                    "#SBATCH -t 8:00:00\n" + \
                    "#SBATCH -A mlodd\n" + \
                    "#SBATCH --mail-type=fail\n" + \
                    f"#SBATCH -e {output_dir}/logs/{label}.err\n" + \
                    f"#SBATCH -o {output_dir}/logs/{label}.out\n" + \
                    "#SBATCH -J " + label + "\n" + \
                    "#SBATCH --export=ALL \n" + \
                    "source /usr/workspace/samuele/spack/share/spack/setup-env.sh \n" + \
                    "spack env activate /usr/workspace/samuele/nexo-spack-env/ \n" + \
                    "source /usr/workspace/samuele/nexo_venv/bin/activate \n" + \
                    "cd " + working_dir + "\n" + \
                    "export STARTTIME=`date +%s`\n" + \
                    "echo Start time $STARTTIME\n" + \
                    "python3 -u " + submit_statement + "\n"\
                    "export STOPTIME=`date +%s`\n" + \
                    "echo Stop time $STOPTIME\n" + \
                    "export DT=`expr $STOPTIME - $STARTTIME`\n" + \
                    "echo CPU time: $DT seconds\n"

        with open(scriptfilename, 'w') as scriptfile:
            scriptfile.write(thescript)

        os.system("sbatch " + scriptfilename)
        