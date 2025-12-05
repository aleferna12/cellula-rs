import argparse
import importlib
import logging
import re
import sys
import os.path


def main():
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(
        prog="scripts",
        description="To see help about each script run 'python -m scripts <script> --help'"
    )
    subparsers = parser.add_subparsers(
        title="scripts",
        required=True,
        description="List of scripts that can be executed via this CLI app"
    )

    # Add root directory to path
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    # Search for the scripts that have added CLI support
    # Include new scripts in this list
    for mod_name in ["fileio",
                     "analyses.make_netgraphs",
                     "analyses.make_trees",
                     "analyses.plot_deaths",
                     "analyses.plot_netgraphs",
                     "analyses.plot_ancestry",
                     "analyses.plot_strategy",
                     "analyses.plot_timeline",
                     "analyses.plot_tree",
                     "analyses.plot_lattice",
                     "analyses.plot_regulation",
                     "analyses.plot_trajectory",
                     "competition.make_competition",
                     "competition.make_templates",
                     "competition.ancestor_template",
                     "colortable.make_colortable",
                     "colortable.print_colortable"]:
        mod = importlib.import_module(f"scripts.{mod_name}")
        subparser = mod.get_parser()
        subparsers.add_parser(
            name=re.search(r"\w+$", mod_name).group(),
            parents=[subparser],
            help=subparser.description,
            add_help=False
        )

    args = parser.parse_args()
    args.run(args)


if __name__ == "__main__":
    main()
