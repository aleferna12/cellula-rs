import argparse
import re
import logging
from pathlib import Path
from functools import partial
from multiprocessing import Pool
from typing import Union
import pandas as pd
import numpy as np
import plotly.express as px
from pathlib import Path
import enlighten

logger = logging.getLogger(__name__)
gene_attrs = ["in_scale_list",
              "reg_threshold_list",
              "reg_w_innode_list",
              "reg_w_regnode_list",
              "out_threshold_list",
              "out_w_regnode_list"]
str_attrs = gene_attrs + ["neighbour_list", "Jneighbour_list"]


def get_parser():
    def run(args):
        celldf = parse_cell_data(args.datadir)
        write_cell_data(celldf, args.outfile)
        logger.info("Finished")

    parser = argparse.ArgumentParser(
        description="Pre-parse cell dataframes into a single CSV file."
    )
    parser.add_argument("datadir",
                        help="Directory containing the cell CSV files")
    parser.add_argument("outfile",
                        help="Output CSV file")
    parser.set_defaults(run=run)
    return parser


def parse_lattice(filepath) -> pd.DataFrame:
    logger.info(f"Parsing dataframe from: '{filepath}'")
    return pd.read_csv(filepath, header=None)


def parse_food_data(datapath,
                    time_filter: Union[list, range] = None,
                    trust_filenames=True,
                    n_processes=1) -> pd.DataFrame:
    df = _parse_dfs(datapath, time_filter, trust_filenames, n_processes=n_processes)
    return df.sort_values("time")


def parse_cell_data(datapath,
                    time_filter: Union[list, range] = None,
                    trust_filenames=True,
                    n_processes=1) -> pd.DataFrame:
    """Parses a CSV file into a dataframe containing cell information.

    The format of this dataframe is essential to API stability, please be mindful of changes to it.
    
    :param n_processes:
    :param time_filter:
    :param datapath:
    :param trust_filenames: Whether to speed up the filtering process by trusting that the file
    names represent the time of the simulation.
    """

    celldf = _parse_dfs(datapath,
                        time_filter,
                        trust_filenames,
                        n_processes,
                        dtype={col: str for col in str_attrs})
    celldf = celldf.sort_values(["time", "sigma"])

    # Can't drop the columns because they might be needed by other functions!
    celldf = celldf.set_index(["time", "sigma"], drop=False)
    celldf.index.names = ["time_i", "sigma_i"]
    return celldf


def parse_grave_data(datapath, time_filter=None, trust_filenames=True, n_processes=1) -> pd.DataFrame:
    gravedf = _parse_dfs(datapath, time_filter, trust_filenames, n_processes=n_processes)
    gravedf = gravedf.sort_values(["time_death", "sigma"])
    return gravedf.set_index(["time_death", "sigma"], drop=False)


def get_time_points(datapath):
    times = set()
    for filepath in Path(datapath).iterdir():
        m = re.search(r"t(\d+).csv", filepath.name)
        if m is not None:
            times.add(int(m.group(1)))
    return sorted(times)


def interspaced_elements(a, n):
    """Gets interspaced elements from a list.

    The elements might not be perfectly interspaced, but 'n' elements are guaranteed to be returned.
    """
    return a[np.round(np.linspace(0, len(a) - 1, n)).astype(int)]


def build_time_filter(time_points, start=0, end=-1, n=None):
    """Selects 'n' interspaced time points in the range ['start', 'end']."""
    if end == -1:
        end = float("inf")

    time_points = np.unique(time_points)
    time_points = time_points[(time_points >= start) & (time_points <= end)]
    if n is None:
        n = len(time_points)
    else:
        n = min(n, len(time_points))
    return interspaced_elements(time_points, n)


def reduce_data(celldf: pd.DataFrame, n):
    times = interspaced_elements(np.unique(celldf["time"]), n)
    logger.info("Reduced cell data frame to contain only the following time-steps:")
    logger.info(times)
    return celldf.loc[celldf["time"].isin(times)]


def write_plot(fig, outputfile):
    logger.info(f"Writing plot to: {outputfile}")
    if outputfile[-5:] == ".html":
        fig.write_html(outputfile)
    else:
        fig.write_image(outputfile)


def write_cell_data(celldf, outfile):
    celldf.to_csv(outfile, index=False)


def write_genomes(celldf, outfile):
    write_cell_data(celldf[["innr", "regnr", "outnr"] + gene_attrs], outfile)


def write_lattice(lattdf, outfile):
    lattdf.to_csv(outfile, header=False, index=False)


def _parse_dfs(datapath, time_filter, trust_filenames, n_processes, **csvkwargs):
    logger.info(f"Parsing dataframe(s) from: '{datapath}'")
    datapath = Path(datapath)

    if datapath.is_file():
        return pd.read_csv(datapath, **csvkwargs)

    filepaths = list(datapath.iterdir())
    if time_filter is not None and trust_filenames:
        filepaths = [fp for fp in filepaths if int(re.search(r"\d+", fp.name).group()) in time_filter]

    manager = enlighten.get_manager(threaded=True, set_scroll=False)
    pbar = manager.counter(total=len(filepaths), desc="CSV files iterated")

    results = []
    pool = Pool(n_processes)
    for fp in filepaths:
        res = pool.apply_async(_read_csv,
                               kwds=dict(filepath=fp, time_filter=time_filter, **csvkwargs),
                               callback=lambda _: pbar.update())
        results.append(res)
    pool.close()
    pool.join()
    manager.stop()
    return pd.concat(r.get() for r in results)


def _read_csv(filepath, time_filter, **csvkwargs):
    df = pd.read_csv(filepath, **csvkwargs)
    if time_filter is None or df.empty or df["time"][0] in time_filter:
        return df
    return None
