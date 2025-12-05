import logging
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)


def restrict_to_range(celldf, par_range):
    old_size = len(celldf)
    celldf = celldf[
        (celldf["grad_conc"].between(par_range.loc["min", "grad_conc"], par_range.loc["max", "grad_conc"])) &
        (celldf["food"].between(par_range.loc["min", "food"], par_range.loc["max", "food"]))
        ]
    logger.info(f"Data frame preserved {len(celldf) / old_size:.1%} of its original size")
    return celldf


def gaussian_kde(celldf):
    datapoints = celldf[["grad_conc", "food"]].values.T
    return stats.gaussian_kde(datapoints)


# TODO: implement or remove out_thres
def make_density_matrix(celldf: pd.DataFrame, grid_size=100, kde=None, n_processes=1):
    # Discard first time step to avoid artificial 0s in chem and food
    celldf = celldf.loc[celldf["time"] != 0]
    # Pre filters very extreme outliers to obtain a better mesh
    # celldf = celldf[
    #     (celldf["grad_conc"].between(*celldf["grad_conc"].quantile([out_thres, 1 - out_thres]))) &
    #     (celldf["food"].between(*celldf["food"].quantile([out_thres, 1 - out_thres])))
    #     ]

    if kde is None:
        kde = gaussian_kde(celldf)

    mesh = np.meshgrid(
        np.linspace(
            celldf["grad_conc"].min(),
            celldf["grad_conc"].max(),
            grid_size
        ),
        np.linspace(
            celldf["food"].min(),
            celldf["food"].max(),
            grid_size
        ),
        indexing="ij"
    )
    mesh = np.reshape(mesh, [2, -1])  # Flattens along one axis and matches cell_counts order

    if n_processes > 1:
        splits = np.array_split(mesh, n_processes, axis=1)
        with Pool(n_processes) as pool:
            densities = pool.map(kde, splits)
        densities = np.concatenate(densities)
    else:
        densities = kde(mesh)
    densities /= np.sum(densities)  # Normalize sum to 1

    chemcats = pd.cut(celldf["grad_conc"], bins=grid_size).cat.categories
    foodcats = pd.cut(celldf["food"], bins=grid_size).cat.categories

    return pd.Series(densities, index=pd.MultiIndex.from_product([chemcats, foodcats]))


def get_parameter_range(celldf):
    return pd.DataFrame(
        {
            "grad_conc": [celldf["grad_conc"].min(), celldf["grad_conc"].max()],
            "food": [celldf["food"].min(), celldf["food"].max()]
        },
        index=["min", "max"]
    )


def filter_kde(celldf, include_thres=0.9, density_matrix=None):
    """ This function filters datapoints in celldf based on kernel density.

    It can be quite slow depending on the size of celldf, and its possible to approximate it by using the
    function observed_parameter_range (that will give you a rectangle where the points are contained).
    """
    og_size = len(celldf)
    celldf = celldf[celldf["time"] != 0]
    if density_matrix is None:
        density_matrix = make_density_matrix(celldf)

    densities = density_matrix.loc[zip(celldf["grad_conc"], celldf["food"])]
    quant = np.quantile(densities, 1 - include_thres)

    filtered = celldf.loc[densities.values > quant]
    logger.info(f"Preserved {len(filtered) / og_size:.1%} of the data points")
    logger.info("Parameter ranges of the filtered data frame are:\n"
                f"\tChem: {celldf['grad_conc'].min()} - {celldf['grad_conc'].max()}\n"
                f"\tFoodp: {celldf['food'].min()} - {celldf['food'].max()}")
    return filtered
