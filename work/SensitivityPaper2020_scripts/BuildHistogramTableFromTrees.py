# Import sys, then tell python where to find the nEXO-specific classes
import sys
import os
sys.path.append('../../modules')

######################################################################
# Check arguments and load inputs
if len(sys.argv) == 5:
        config_file = sys.argv[1]
        label = sys.argv[2]
        path_to_trees = sys.argv[3]
        output_dir = sys.argv[4]
        if not os.path.exists(output_dir):
                sys.exit('\nERROR: path to output_dir does not exist\n')
else:
        print('\n\nERROR: BuildHistogramsFromTrees.py requires 3 arguments')
        print('Usage:')
        print('\tpython BuildHistogramsFromTrees.py ' + \
                '<config_file> <label> <path/to/merged/ROOT/trees> <output_dir>')
        sys.exit('\n')
######################################################################

# Import the nEXO sensitivity classes
import nEXOFitWorkspace

workspace = nEXOFitWorkspace.nEXOFitWorkspace( config = config_file )

output_hdf5_filename = '{}/Baseline2019_Histograms_{}.h5'.format(output_dir,label)

workspace.CreateHistogramsFromRawTrees( path_to_trees = path_to_trees, \
                             output_hdf5_filename = output_hdf5_filename )


