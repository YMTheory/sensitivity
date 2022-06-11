# Import sys, then tell python where to find the nEXO-specific classes
import sys
import argparse
import pathlib
sys.path.append('../../modules')


def get_parser():
    parser = argparse.ArgumentParser(description='Build HDF5 Histogram Table from Raw TTrees')
    parser.add_argument("config_file", type=pathlib.Path, help='nEXOFitWorkspace configuration file')
    parser.add_argument("label", type=str, help='Label')
    parser.add_argument("path_to_trees", type=pathlib.Path, help='path/to/merged/ROOT/trees')
    parser.add_argument("output_dir", type=pathlib.Path, help='Path where the HDF5 file will be saved')
    parser.add_argument("-r", "--resolution_factor", type=float, default=0, help='Energy resolution smearing factor')
    parser.add_argument("-d", "--dnn_smoothing_factor", type=float, default=0, help='DNN smoothing factor')
    return parser


if __name__ == "__main__":
    arg_parser = get_parser()
    args = arg_parser.parse_args()
    for arg in vars(args):
        print(arg, getattr(args, arg))

    # Import the nEXO sensitivity classes
    import nEXOFitWorkspace

    workspace = nEXOFitWorkspace.nEXOFitWorkspace(config=args.config_file)

    output_hist_basename = args.output_dir / pathlib.Path(f'Baseline2019_Histograms_{args.label}')

    workspace.CreateHistogramsFromRawTrees(path_to_trees=str(args.path_to_trees),
                                           output_hist_basename=output_hist_basename,
                                           resolution_factor=args.resolution_factor,
                                           dnn_smoothing_factor=args.dnn_smoothing_factor
                                           )


