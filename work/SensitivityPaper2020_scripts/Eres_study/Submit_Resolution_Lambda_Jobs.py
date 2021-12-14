#!/usr/local/bin/python

import os

caldate = '21_03_11'
compdate = '21_03_01'

execdir = "/p/lustre2/czyz1/nexo_sensitivity/work/SensitivityPaper2020_scripts"
outputdir = "/p/lustre2/nexouser/czyz1/output"
config_loc = "/p/lustre2/czyz1/nexo_sensitivity/work/config/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1.yaml"
USE_WILKS = 'True'

bkg_shape_err = 0.00000001
num_datasets = 100
num_hypotheses = 15

base = "Sensitivity_run_trial_"

for database in ['024']:
# for database in ['023', '024', '025']:
    date = caldate + '_DNN1_' + database
    lamdate = '21_03_01_DNN1_' + database
    comp_loc = "/p/lustre2/nexouser/czyz1/workdir/components_tables/{}/ComponentsTable_D-{}_Energy_Res=".format(compdate, database)
    # comp_loc = '/p/vast1/nexo/sensitivity2020/pdfs/component_tables/ComponentsTable_D-025_merged-v10b_Optimized_DNN_Standoff_Binning_version1.h5'
    crit_lam_loc = "/p/lustre2/nexouser/czyz1/workdir/lambda/{}".format(lamdate)

    if database == '024':
        # resolution_list = ['0.008', '0.009', '0.01', '0.011', '0.012', '0.013', '0.014', '0.015', '0.016', '0.017', '0.018']
        # resolution_list = ['0.008', '0.01', '0.012', '0.014', '0.016', '0.018']
        resolution_list = ['0.009', '0.011', '0.013', '0.015', '0.017']
    else:
        resolution_list = ['0.008']

    for resolution in resolution_list:
        for num in range(0, 100):
            basename = base + str(num)
            if not os.path.exists(outputdir + "/Sub/" + date):
                    os.makedirs(outputdir + "/Sub/" + date)
            scriptfilename = outputdir + "/Sub/" + date + "/" + base + "_" + resolution + "_" + str(num) + ".sub"
            os.system("rm -f " + scriptfilename)
            if not os.path.exists(outputdir + "/Out/" + date):
                    os.makedirs(outputdir + "/Out/" + date)
            outfilename = outputdir + "/Out/" + date + "/" + base + "_" + resolution + "_" + str(num) + ".out"
            os.system("rm -f " + outfilename)
            errfilename = outputdir + "/Err/" + base  + "_" + resolution + "_" + str(num) + ".err"
            os.system("rm -f " + errfilename)

            thescript = "#!/bin/bash\n" + \
                    "#SBATCH -t 96:00:00\n" + \
                    "#SBATCH -A nuphys\n" + \
                    "#SBATCH -e " + outfilename + "\n" + \
                    "#SBATCH -o " + outfilename + "\n" + \
                    "#SBATCH --mail-type=fail\n" + \
                    "#SBATCH -J " + base + "\n" + \
                    "#SBATCH --export=ALL \n" + \
                    "#SBATCH -p pbatch\n" + \
                    "source /usr/gapps/nexo/setup.sh \n" + \
                    "source /g/g12/czyz1/nexoenv/bin/activate \n" + \
                    "cd " + execdir + "\n" + \
                    "export STARTTIME=`date +%s`\n" + \
                    "echo Start time $STARTTIME\n" + \
                    "python3 " + execdir + "/Eres_study/Resolution_Lambda_Sensitivity.py {} {} {} {} {} {} {} {} {} {} {}" \
                                           "\n".format(num, bkg_shape_err, num_datasets, outputdir, num_hypotheses,
                                               date, config_loc, comp_loc, crit_lam_loc, resolution, USE_WILKS) + \
                    "export STOPTIME=`date +%s`\n" + \
                    "echo Stop time $STOPTIME\n" + \
                    "export DT=`expr $STOPTIME - $STARTTIME`\n" + \
                    "echo CPU time: $DT seconds\n"

            scriptfile = open(scriptfilename, 'w')
            scriptfile.write(thescript)
            scriptfile.close()

            os.system("sbatch " + scriptfilename)
