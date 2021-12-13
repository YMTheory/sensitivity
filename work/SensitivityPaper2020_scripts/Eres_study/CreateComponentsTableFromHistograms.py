# Import syn, then tell python where to find the nEXO-specific classes
import sys
import os
sys.path.append('../../../modules')

geometry_tags = ['D-023','D-024','D-025']
#geometry_tags = ['D-025']
date = '21_03_01'
for geometry_tag in geometry_tags:
    for resolution in ['0.008', '0.009', '0.01', '0.011', '0.012', '0.013', '0.014', '0.015', '0.016', '0.017', '0.018']:
    # for resolution in ['0.008']:
        input_histogram_file = '/p/lustre2/nexouser/czyz1/workdir/histogram_files/' + date + '/Baseline2019_Histograms_Energy_Res={' \
                               '}.h5'.format(resolution)
        label = 'Energy_Res={}'.format(resolution)
        config_file = '../../config/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1.yaml'
        output_dir = '/p/lustre2/nexouser/czyz1/workdir/components_tables/' + date + '/'
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        import nEXOFitWorkspace

        workspace = nEXOFitWorkspace.nEXOFitWorkspace(config=config_file)

        workspace.CreateComponentsTableFromMaterialsDB( geometry_tag=geometry_tag,\
                                        histograms_file = input_histogram_file,\
                                        label=label, output_dir=output_dir)

