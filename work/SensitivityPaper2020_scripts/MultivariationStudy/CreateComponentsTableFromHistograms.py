# Import syn, then tell python where to find the nEXO-specific classes
import sys
import os

sys.path.append('../../../modules')

geometry_tags = ['D-024']
date = '21_11_24'
dnn_factors = [0., 0.15, 0.2]

for geometry_tag in geometry_tags:
    for dnn_factor in dnn_factors:
        # for resolution in ['0.008']:
        label = f'DNN_factor={dnn_factor}'
        input_histogram_file = '/Users/sangiorgio1/OneDrive - LLNL/nEXO/Simulations/sensitivity/sensitivity2020/multi_variations_study/' + date + \
                               f'/Baseline2019_Histograms_{label}.h5'
        config_file = '/Users/sangiorgio1/OneDrive - LLNL/nEXO/Simulations/sensitivity/sensitivity2020/multi_variations_study/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1.yaml'
        output_dir = '/Users/sangiorgio1/scratch/' + date + '/'
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        import nEXOFitWorkspace

        workspace = nEXOFitWorkspace.nEXOFitWorkspace(config=config_file)

        workspace.CreateComponentsTableFromMaterialsDB(geometry_tag=geometry_tag,
                                                       histograms_file=input_histogram_file,
                                                       label=label, output_dir=output_dir)
