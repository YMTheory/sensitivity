# Import sys, then tell python where to find the nEXO-specific classes
import sys
import os
sys.path.append('../../modules')

######################################################################
# Check arguments and load inputs
if len(sys.argv) == 6:
        config_file = sys.argv[1]
        label = sys.argv[2]
        input_histogram_file = sys.argv[3]
        output_dir = sys.argv[4]
        tag = sys.argv[5]
        if not os.path.exists(output_dir):
                sys.exit('\nERROR: path to output_dir does not exist\n')
else:
        print('\n\nERROR: CreateComponentsTableFromHistograms.py requires 5 arguments')
        print('Usage:')
        print('\tpython CreateComponentsTableFromHistograms.py ' + \
                '<config_file> <label> <histograms_file> <output_dir> <geometry_tag>')
        sys.exit('\n')
######################################################################


import nEXOFitWorkspace

workspace = nEXOFitWorkspace.nEXOFitWorkspace(config=config_file)

workspace.CreateComponentsTableFromMaterialsDB(geometry_tag=tag,\
                               histograms_file=input_histogram_file,\
                               label=label, output_dir=output_dir)
