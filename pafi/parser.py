import argparse
import pafi.plots as plots

def parse_cli_args():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(help='', required=True, title='Available modes',)

    ## NEB MODE
    parser_NEB = subparsers.add_parser('neb', aliases=["NEB"], help='Plot the energy profile from a NEB log.lammps file.')
    parser_NEB.add_argument('log_lammps', 
                        help="log.lammps file from a NEB run")
    parser_NEB.add_argument('-f', '--final', action="store_true",
                        help="Only display the final stage of the relaxation. Otherwise, profiles at the initial and pre-climb stages are also shown.")
    parser_NEB.add_argument('-C', '--cut-at-image', default=-1, metavar="N",
                        help="If N is specified, only the first N images will be shown.")
    parser_NEB.add_argument('-n', '--normalize-r', action="store_true",
                        help="Normalize r such that its maximum value is 1.0")
    parser_NEB.add_argument('-N', '--image-number', action="store_true",
                        help="Show image number instead of reaction coordinate on the x axis")
    parser_NEB.add_argument('-E', '--energy_diff', action="store_true", 
                        help="Add a subplot showing the spacing of points on the Energy axis. Useful to spot large gaps in energy that would cause a large error for PAFI integration.")
    parser_NEB.add_argument('-R', '--rcdiff', action="store_true", 
                        help="Add a subplot showing the spacing of points on the r axis. Useful to spot large gaps in reaction coordinate.")
    parser_NEB.add_argument('--print-profile', action="store_true", 
                        help="Print NEB profile raw data to terminal.")
    parser_NEB.set_defaults(func=plots.plot_neb_results)

    ## PAFI MODE
    parser_pafi = subparsers.add_parser('pafi', aliases=["PAFI"], help='Plot the free energy profile from a PAFI run using a raw_ensemble_output file.')
    parser_pafi.add_argument('file', 
                        help="A PAFI raw_ensemble_output file, or a glob expression that matches a list of raw_ensemble_output files.")
    parser_pafi.add_argument('--mode', choices="sm",
                        help="either plot a single or multiple raw files.")
    parser_pafi.add_argument('--H0', '--zero', action="store",
                        help="A value for 0K, e.g obtained from NEB.")
    parser_pafi.add_argument('--no-fit', action="store_false",
                        help="Disable fitting DeltaS on the harmonic domain.")
    parser_pafi.add_argument('--hist', action="store_true",
                        help="Plot an histogram of dF.")
    parser_pafi.add_argument('--std', action="store_true",
                        help="Plot std for each hyperplane.")
    parser_pafi.add_argument('--harmonic-limit', default=100, type=float,
                        help="The upper temperature limit above which the entropy is considered non-linear. System dependant.")
    parser_pafi.set_defaults(func=plots.plot_pafi_results)

    # COMMON TO ALL MODES 
    for p in (parser_pafi, parser_NEB): 
        p.add_argument('-t', '--terminal', action="store_true", 
                            help="Experimental: use plotext as backend for display in terminal.")
        p.add_argument('-s', '--save', 
                        help='name of a file destination (e.g. "plot.pdf", "plot.png")')
        p.add_argument('-q', '--quiet', action="store_true",
                        help="Do not open matplotlib popup nor display in terminal.")
        p.add_argument('--size', default=4, 
                        help="Scaling factor of the figure passed to matplotlib (default is 4).")

    return parser.parse_args()