import pafi.parser
import pafi.plots
from pafi import pafi

def cli():
    args = pafi.parser.parse_cli_args()
    args.func(args)
