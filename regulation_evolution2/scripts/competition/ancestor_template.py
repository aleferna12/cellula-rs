import argparse
import logging
import pandas as pd
from scripts.fileio import *
from scripts.analyses.make_trees import make_trees
from scripts.competition.make_templates import template_from_cells

logger = logging.getLogger(__name__)


def get_parser():
    def run(args):
        celldf = parse_cell_data(args.datadir)
        ancestor_template(celldf=celldf,
                          outfile=args.outfile,
                          food=args.food)
        logger.info("Finished")

    parser = argparse.ArgumentParser(
        description="Retrieve the ancestor of a simulation and make it a template."
    )
    parser.add_argument(
        "datadir",
        help="Directory containing the cell data frames from which the ancestor will be extracted"
    )
    parser.add_argument("outfile", help="Output file for the ancestor template")
    parser.add_argument("-f", "--food", help="How much food does each cell start with")
    parser.set_defaults(run=run)
    return parser


def retrieve_ancestor(celldf: pd.DataFrame):
    trees = make_trees(celldf, stop_mrca=True)
    if len(trees) != 1:
        raise ValueError("population doesn't have a single ancestor")

    ancestor_node = trees[0]
    return celldf.loc[(ancestor_node.time, int(ancestor_node.name))]


def ancestor_template(celldf: pd.DataFrame, outfile, food=-1):
    logger.info("Retrieving ancestor")
    adf = retrieve_ancestor(celldf).to_frame().T
    logger.info(f"Ancestor had sigma {adf.iloc[0]['sigma']} and lived at MCS {adf.iloc[0]['time']}")
    templatedf = template_from_cells(adf, food=food)
    write_cell_data(templatedf, outfile)
    return templatedf
