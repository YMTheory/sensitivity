# Import sys, then tell python where to find the nEXO-specific classes
import sys
import os
import glob

sys.path.append('../../../modules')

geometry_tags = ['D-024']
date = '22_03_30_fineBinning'
dnn_factors = [0., 0.15, 0.177, 0.2, 0.25]
# dnn_factors = [0.15, ]
resolutions = [0.008, 0.01, 0.011, 0.012, 0.014]
# resolutions = [0.011]

for geometry_tag in geometry_tags:
    for resolution in resolutions:
        for dnn_factor in dnn_factors:
            label = f'DNN_factor={dnn_factor}_ERes={resolution}_fineBinning'
            input_histogram_file = '/p/lustre2/nexouser/samuele/multivarstudy/' + date + \
                                   f'/Baseline2019_Histograms_{label}.pkl.gz'
            config_file = '/g/g92/samuele/nEXO/sensitivity/work/SensitivityPaper2020_scripts/MultivariationStudy/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1_fineBinning.yaml'
            output_dir = '/p/lustre2/nexouser/samuele/multivarstudy/' + date + '/'
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
                
            if os.path.exists(f"{output_dir}ComponentsTable_{geometry_tag}_{label}.pkl.gz"):
                print(f"Output for {label} already exists. Skipping.")
                continue
                             
            print(f"Processing {label}...")
            
            import nEXOFitWorkspace

            workspace = nEXOFitWorkspace.nEXOFitWorkspace(config=config_file)

            workspace.CreateComponentsTableFromMaterialsDB(geometry_tag=geometry_tag,
                                                           histograms_file=input_histogram_file,
                                                           label=label, output_dir=output_dir)
