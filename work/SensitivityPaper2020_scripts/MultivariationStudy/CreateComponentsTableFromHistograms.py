# Import sys, then tell python where to find the nEXO-specific classes
import sys
import os

sys.path.append('../../../modules')

geometry_tags = ['D-024']
date = '22_01_14'
# dnn_factors = [0., 0.15, 0.2, 0.25]
dnn_factors = [0.177, ]
resolutions = [0.008, 0.01, 0.012, 0.014]
#resolutions = [0.014]

for geometry_tag in geometry_tags:
    for resolution in resolutions:
        for dnn_factor in dnn_factors:
            label = f'DNN_factor={dnn_factor}_ERes={resolution}'
            input_histogram_file = '/p/lustre2/nexouser/samuele/multivarstudy/' + date + \
                                   f'/Baseline2019_Histograms_{label}.h5'
            config_file = '/g/g92/samuele/nEXO/sensitivity/work/SensitivityPaper2020_scripts/MultivariationStudy/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1.yaml'
            output_dir = '/p/lustre2/nexouser/samuele/multivarstudy/' + date + '/'
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            import nEXOFitWorkspace

            workspace = nEXOFitWorkspace.nEXOFitWorkspace(config=config_file)

            workspace.CreateComponentsTableFromMaterialsDB(geometry_tag=geometry_tag,
                                                           histograms_file=input_histogram_file,
                                                           label=label, output_dir=output_dir)
