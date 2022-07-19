# Import sys, then tell python where to find the nEXO-specific classes
import sys
import argparse
import pathlib
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
    parser.add_argument('-p', '--path_to_data', type=pathlib.Path, \
                        default = None, help = 'path to the merged ROOT data')
    parser.add_argument('-o', '--output_dir', type=pathlib.Path, \
                        default = None, help = 'output directory')
    parser.add_argument('-r', '--resolution_factor', type=float, \
                        default = 0, help = 'smearing factor for energy resolution')
    parser.add_argument('-d', '--dnn_smoothing_factor', type=float, \
                        default=0, help='DNN smoothing factor')
    return parser
######################################################################


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    for arg in vars(args):
        print(arg, getattr(args, arg))
     
    for arg, val in vars(args).items():
        if val is None:
           parser.print_help()
           sys.exit()   

    workspace = nEXOFitWorkspace.nEXOFitWorkspace( config = args.config )


    output_hist_basename = args.output_dir / pathlib.Path(f'Baseline2019_Histograms_{args.label}')

    workspace.CreateHistogramsFromRawTrees(path_to_trees=str(args.path_to_trees),
                                           output_hist_basename=output_hist_basename,
                                           resolution_factor=args.resolution_factor,
                                           dnn_smoothing_factor=args.dnn_smoothing_factor
                                           )


