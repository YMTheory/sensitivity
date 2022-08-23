# Import sys, then tell python where to find the nEXO-specific classes
import sys
import os
import argparse
sys.path.append('../../modules')
import nEXOFitWorkspace 

######################################################################
# Parse arguments
def get_parser():
    parser = argparse.ArgumentParser( description = \
                                      'Create components table from histlite histograms')
    parser.add_argument('-c', '--config', type=str, default = None, help = 'config file')
    parser.add_argument('-l', '--label', type=str, default = None, help = 'user-defined label')
    parser.add_argument('-f', '--histfile', type=str, default = None, help = 'file containing histograms')
    parser.add_argument('-o', '--output_dir', type=str, default = None, help = 'output directory')
    parser.add_argument('-t', '--tag', type=str, default = None, help = 'geometry tag')
    return parser
######################################################################


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    print(args)

    for arg, val in vars(args).items():
        if val is None:
           parser.print_help()
           sys.exit()   
 

    workspace = nEXOFitWorkspace.nEXOFitWorkspace( config = args.config)
    
    workspace.CreateComponentsTableFromMaterialsDB( geometry_tag = args.tag,\
                                   histograms_file = args.histfile,\
                                   label = args.label, output_dir = args.output_dir)
