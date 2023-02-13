#!/usr/local/bin/python
import hashlib
import os

execdir = "/g/g92/samuele/nEXO/sensitivity/work/SensitivityPaper2020_scripts/SiteSelectionStudy"
outputdir = "/p/lustre2/nexouser/samuele/siteselection"
executable_name = 'Compute90PercentLimit_WilksApprox.py'
components_table_dir = '/p/lustre2/nexouser/samuele/multivarstudy/ComponentsTables'
components_table_basename = 'ComponentsTable_D-024'
config_file = '/g/g92/samuele/nEXO/sensitivity/work/SensitivityPaper2020_scripts/SiteSelectionStudy/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1_Xe137.yaml'

# one full calculation for 5000 toys takes about 4-5 hours
bkg_shape_err = 0.
num_datasets = 6000
num_jobs = 100
jobs_offset = 0
num_datasets_per_job = int(num_datasets / num_jobs)

label = "LNGS_RV"
base = "Run_siteselection_"

components_table = f'{components_table_basename}_DNN_factor=0.0_ERes=0.008.h5'

for iter_num in range(jobs_offset, num_jobs + jobs_offset):

    basename = base + f"{iter_num}_" + label
    scriptfilename = f'{outputdir}/jobs/{basename}.sh'

    os.makedirs(f'{outputdir}/logs/', exist_ok=True)
    os.makedirs(f'{outputdir}/jobs/', exist_ok=True)
    os.makedirs(f'{outputdir}/output/', exist_ok=True)
    os.system(f"rm -f {scriptfilename}")
    os.system(f"rm -f {outputdir}/logs/{basename}.*")

    submit_statement = f'{executable_name} {iter_num} {components_table_dir}/{components_table} ' + \
                       f'{outputdir}/output/ -c {config_file} -n {num_datasets_per_job} -l {label}' 

    thescript = "#!/bin/bash\n" + \
                "#SBATCH -t 12:00:00\n" + \
                "#SBATCH -A mlodd\n" + \
                f"#SBATCH -e {outputdir}/logs/{basename}.err\n" + \
                f"#SBATCH -o {outputdir}/logs/{basename}.out\n" + \
                "#SBATCH --mail-type=fail\n" + \
                f"#SBATCH -J {label}\n" + \
                "#SBATCH --export=ALL \n" + \
                "source /usr/workspace/samuele/nexo_venv/bin/activate \n" + \
                "cd " + execdir + "\n" + \
                "export STARTTIME=`date +%s`\n" + \
                "echo Start time $STARTTIME\n" + \
                f"python3 -u {submit_statement}\n" + \
                "export STOPTIME=`date +%s`\n" + \
                "echo Stop time $STOPTIME\n" + \
                "export DT=`expr $STOPTIME - $STARTTIME`\n" + \
                "echo CPU time: $DT seconds\n"

    with open(scriptfilename, 'w') as scriptfile:
        scriptfile.write(thescript)

    os.system("sbatch " + scriptfilename)
