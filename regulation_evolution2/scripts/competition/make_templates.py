import argparse
import logging
import numpy as np
import pandas as pd
from typing import List
from ..fileio import *

logger = logging.getLogger(__name__)


def get_parser():
    def run(args):
        celldfs = [parse_cell_data(path) for path in args.inputfiles]
        sample_template(celldfs, args.outfile, args.number, args.food)
        logger.info("Finished")

    parser = argparse.ArgumentParser(
        description="Create cell templates used to generate competition files with "
                    "'make_competition'."
    )
    parser.add_argument("outfile",
                        help="CSV output template file")
    parser.add_argument("-i",
                        "--inputfiles",
                        nargs='+',
                        required=True,
                        help="Space-delimited list of cell dataframe CSV files extracted from simulations. "
                             "One or more cells from each of these files will be added to the "
                             "output template file")
    parser.add_argument("-n",
                        "--number",
                        default=1,
                        type=int,
                        help="Number of cells to sample from each CSV file (default: %(default)s)")
    parser.add_argument("-f",
                        "--food",
                        default=-1,
                        type=int,
                        help="The amount of food the will be reassigned to each cell in the "
                             "final template. Use -1 to keep the food the cell had (default: %(default)s)")
    parser.set_defaults(run=run)
    return parser


def template_from_cells(celldf, food=-1):
    """Resets attributes to turn a data frame of cells into a template."""
    celldf = celldf.copy()
    celldf["time"] = 0
    celldf["time_since_birth"] = 0
    celldf["last_meal"] = 0
    celldf["times_divided"] = 0
    celldf["dividecounter"] = 0
    celldf["tau"] = 1
    if food != -1:
        celldf["food"] = food
    return celldf


def sample_template(celldfs: List[pd.DataFrame], outputfile, n=10, food=200):
    """Sample cells from cell data frames and returns a template."""
    templates = []
    for i, celldf in enumerate(celldfs):
        celldf = celldf.sample(n)
        celldf.group = i
        templates.append(celldf)
    tdf = template_from_cells(pd.concat(templates), food)
    tdf.to_csv(outputfile, index=False)


