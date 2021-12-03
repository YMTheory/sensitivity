#!/usr/local/bin/python
import hashlib
import os
import itertools

execdir = "/g/g92/samuele/nEXO/sensitivity/work/SensitivityPaper2020_scripts/MultivariationStudy/"
outputdir = "/p/lustre2/nexouser/samuele/multivarstudy/"
executable_name = 'Compute90PercentLimit_WilksApprox_MultivarStudy.py'
components_table_dir = '/p/lustre2/nexouser/samuele/multivarstudy/ComponentsTables/'
components_table_basename = 'ComponentsTable_D-024_'
config_file = '/g/g92/samuele/nEXO/sensitivity/work/SensitivityPaper2020_scripts/MultivariationStudy/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1.yaml'

dnn_factors = [0., 0.15, 0.2]
# xe137_scale_factors = [1.,0.01,0.1,0.3,3.,10.,30.,100.]
xe137_scale_factors = [1., ]
rn222_scale_factors = [1., ]
bkg_scale_factors = [1., ]
energy_res_factors = [0.008, ]

bkg_shape_err = 0.
num_datasets = 5000
num_jobs = 100
jobs_offset = 0
num_datasets_per_job = int(num_datasets / num_jobs)

base = "Run_multivar_"

for dnn_scale_factor, xe137_scale_factor, rn222_scale_factor, bkg_scale_factor, energy_res \
        in itertools.product(dnn_factors, xe137_scale_factors, rn222_scale_factors, bkg_scale_factors,
                             energy_res_factors):
    s = f'Xe137:{xe137_scale_factor:0>4.4f} ' + \
        f'Rn222:{rn222_scale_factor:0>4.4f} ' + \
        f'DNN:{dnn_scale_factor:0>4.4f} ' + \
        f'Bkg:{bkg_scale_factor:0>4.4f} ' + \
        f'ERes:{energy_res:0>4.4f}'
    s_hash = hashlib.md5(s.encode('utf-8')).hexdigest()[:6].upper()

    components_table = f'{components_table_basename}_DNN_factor={dnn_scale_factor}_ERes={energy_res}.h5'

    print(f'\n\nSubmitting jobs for hash {s_hash} with parameters: {s}.\n\n')
    for iter_num in range(jobs_offset, num_jobs + jobs_offset):

        scriptfilename = outputdir + base + s_hash + f"_{iter_num:03}.sh"
        os.system("rm -f " + scriptfilename)
        outfile_basename = outputdir + base + s_hash + f"_{iter_num:03}"
        os.system("rm -f " + outfilename)

        submit_statement = executable_name + f'{iter_num} {components_table_dir}/{components_table} ' + \
                           f'{outputdir} -c {config_file} -n {num_datasets_per_job} -e {energy_res} ' + \
                           f'-d {dnn_scale_factor} -b {bkg_scale_factor} -x {xe137_scale_factor} -r {rn222_scale_factor}'

        thescript = "#!/bin/bash\n" + \
                    "#SBATCH -t 14:00:00\n" + \
                    "#SBATCH -A nuphys\n" + \
                    f"#SBATCH -e {outfile_basename}.err\n" + \
                    f"#SBATCH -o {outfile_basename}.out\n" + \
                    "#SBATCH --mail-type=fail\n" + \
                    f"#SBATCH -J {s_hash}\n" + \
                    "#SBATCH --export=ALL \n" + \
                    "source /usr/gapps/nexo/setup.sh \n" + \
                    "source /usr/workspace/samuele/nexo_venv/bin/activate \n" + \
                    "cd " + execdir + "\n" + \
                    "export STARTTIME=`date +%s`\n" + \
                    "echo Start time $STARTTIME\n" + \
                    "python3" + submit_statement + \
                    "export STOPTIME=`date +%s`\n" + \
                    "echo Stop time $STOPTIME\n" + \
                    "export DT=`expr $STOPTIME - $STARTTIME`\n" + \
                    "echo CPU time: $DT seconds\n"

        scriptfile = open(scriptfilename, 'w')
        scriptfile.write(thescript)
        scriptfile.close()

        os.system("sbatch " + scriptfilename)
