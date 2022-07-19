# Import sys, then tell python where to find the nEXO-specific classes
import sys
import os
import argparse
sys.path.append('../../modules')

# Import the nEXO sensitivity classes
import nEXOFitWorkspace


######################################################################
# Parse arguments
def get_parser():
    parser = argparse.ArgumentParser( description = \
                                      'Build the table of histlite histograms' +
                                      ' from the processed & merged ROOT data')
    parser.add_argument('-c', '--config', type=str, \
                        default = None, help = 'config file')
    parser.add_argument('-l', '--label', type=str, \
                        default = None, help = 'user-defined label')
    parser.add_argument('-p', '--path_to_data', type=str, \
                        default = None, help = 'path to the merged ROOT data')
    parser.add_argument('-o', '--output_dir', type=str, \
                        default = None, help = 'output directory')
    parser.add_argument('-r', '--resolution_factor', type=str, \
                        default = None, help = 'smearing factor for energy resolution')
    return parser
######################################################################


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    print(args)
     
    for arg, val in vars(args).items():
        if val is None and not arg is "resolution_factor":
           parser.print_help()
           sys.exit()   

    workspace = nEXOFitWorkspace.nEXOFitWorkspace( config = args.config )

    output_hdf5_filename = '{}/Histograms_{}.h5'.format( args.output_dir, args.label)

    workspace.CreateHistogramsFromRawTrees( path_to_trees = args.path_to_data, \
                             output_hdf5_filename = output_hdf5_filename, \
                             resolution_factor = args.resolution_factor )


