#!/usr/local/bin/python

import os
execdir = "/p/lustre2/czyz1/nexo_sensitivity/work/SensitivityPaper2020_scripts"
outputdir = "/p/lustre2/nexouser/czyz1/output"
config_loc = "/p/lustre2/czyz1/nexo_sensitivity/work/config/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1_v9wAr42.yaml"
comp_loc = "/p/lustre2/nexouser/czyz1/workdir/components_tables/ComponentsTable_D-024_wAr42_Energy_Res="

bkg_shape_err = 0.00000001
num_datasets = 10
num_hypotheses = 15

base = "Sensitivity_run_test_trial_"
date = "20_11_30_DNN1_024"

for num in range(0, 250):

        basename = base + str(num)
        if not os.path.exists(outputdir + "/Sub/" + date):
                os.makedirs(outputdir + "/Sub/" + date)
        scriptfilename = outputdir + "/Sub/"+ date + "/" + base + str(num) + ".sub"
        os.system("rm -f " + scriptfilename)
        if not os.path.exists(outputdir + "/Out/"+ date):
                os.makedirs(outputdir + "/Out/"+ date)
        outfilename = outputdir + "/Out/" + date + "/" + base + str(num) + ".out"
        os.system("rm -f " + outfilename)
        errfilename = outputdir + "/Err/" + base + str(num) + ".err"
        os.system("rm -f " + errfilename)

        thescript = "#!/bin/bash\n" + \
                "#SBATCH -t 48:00:00\n" + \
                "#SBATCH -A nuphys\n" + \
                "#SBATCH -e " + outfilename + "\n" + \
                "#SBATCH -o " + outfilename + "\n" + \
                "#SBATCH --mail-type=fail\n" + \
                "#SBATCH -J " + base + "\n" + \
                "#SBATCH --export=ALL \n" + \
                "source /usr/gapps/nexo/setup.sh \n" + \
                "source /g/g12/czyz1/nexoenv/bin/activate \n" + \
                "cd " + execdir + "\n" + \
                "export STARTTIME=`date +%s`\n" + \
                "echo Start time $STARTTIME\n" + \
                "python3 " + execdir + "/Eres_study/Resolution_Wilks_Sensitivity.py {} {} {} {} {} {} {} {}\n".format(\
                               num, bkg_shape_err, num_datasets, outputdir, num_hypotheses, date, config_loc, comp_loc) + \
                "export STOPTIME=`date +%s`\n" + \
                "echo Stop time $STOPTIME\n" + \
                "export DT=`expr $STOPTIME - $STARTTIME`\n" + \
                "echo CPU time: $DT seconds\n"

        scriptfile = open(scriptfilename, 'w')
        scriptfile.write(thescript)
        scriptfile.close()

        os.system("sbatch " + scriptfilename)

