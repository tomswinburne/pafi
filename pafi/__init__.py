import pafi.parser
import pafi.plots

def cli():
    args = pafi.parser.parse_cli_args()
    args.func(args)
