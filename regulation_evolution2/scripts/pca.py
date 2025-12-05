import sys
import numpy as np
import pandas as pd
import plotly.express as px
import colorir as cl
from pathlib import Path
from logging import getLogger
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scripts.calculate_adh import *
from scripts.fileio import *
from scripts.sweep import sweep_cells
from scripts.analyses.make_netgraphs import observed_parameter_range, restrict_to_range
from scripts.analyses.plot_ancestry import get_mrcass

logger = getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO)
    celldatadirs = Path(sys.argv[1]).expanduser().iterdir()
    celldfs = {}
    par_range = None
    for datapath in celldatadirs:
        runname = datapath.name
        celldir = datapath / "celldata_"
        tf = build_time_filter(get_time_points(celldir), 10000, n=100)
        celldfs[runname] = parse_cell_data(celldir, time_filter=tf, n_processes=15)
        new_par_range, _ = observed_parameter_range(
            celldfs[runname],
            include_thres=0.5
        )
        if par_range is None:
            par_range = new_par_range
        else:
            par_range.loc["min", "grad_conc"] = min(par_range.loc["min", "grad_conc"],
                                                    new_par_range.loc["min", "grad_conc"])
            par_range.loc["max", "grad_conc"] = max(par_range.loc["max", "grad_conc"],
                                                    new_par_range.loc["max", "grad_conc"])
            par_range.loc["min", "food"] = min(par_range.loc["min", "food"],
                                               new_par_range.loc["min", "food"])
            par_range.loc["max", "food"] = max(par_range.loc["max", "food"],
                                               new_par_range.loc["max", "food"])

        if len(sys.argv) >= 6:
            times = interspaced_elements(np.unique(celldfs[runname]["time"]), int(sys.argv[5]))
            celldfs[runname] = celldfs[runname][celldfs[runname]["time"].isin(times)]

    simdf = concatenate_sims(celldfs)
    simdf = restrict_to_range(simdf, par_range)
    logger.info(f"Chem range before aggregation is: {simdf['grad_conc'].min()} - {simdf['grad_conc'].max()}")
    logger.info(f"Food range before aggregation is: {simdf['food'].min()} - {simdf['food'].max()}")
    datadf = make_datadf(simdf, int(sys.argv[2]), int(sys.argv[3]))

    parspace = make_parspacedf(datadf.dropna())
    logger.info("This is the number of simulations that contain data for each parameter coordinate:")
    logger.info(parspace.to_string())

    reprcelldf = concatenate_sims({k: celldf.iloc[0].to_frame().T for k, celldf in celldfs.items()}).set_index("sim")
    chemindex = datadf.index.levels[1]
    foodindex = datadf.index.levels[2]
    sweepdf = sweep_cells(reprcelldf,
                          min_chem=chemindex[0].right,
                          max_chem=chemindex[-1].right,
                          step_chem=chemindex[-1].length,
                          min_foodp=foodindex[0].right / 204,
                          max_foodp=foodindex[-1].right / 204,
                          step_foodp=foodindex[-1].length / 204,
                          mcss=50,
                          reset=True)
    datadf = handle_missing_data(datadf, sweepdf)
    pcadf = make_pcadf(datadf, 24, 12, range(1, 9))

    results, pca = run_pca(pcadf)
    logger.info(f"Cumulative explained variance ratio is: {pca.explained_variance_ratio_.cumsum()}")
    n_sims = datadf.index.get_level_values('sim').nunique()
    colormap = {str(i): p for i, p in enumerate(cl.StackPalette.load("carnival").resize(n_sims))}
    write_plot(px.scatter(
        results,
        x="pc1",
        y="pc2",
        color=results.index,
        color_discrete_map=colormap
    ), sys.argv[4])
    pass


def handle_missing_data(datadf, sweepdf):
    fill_datadf(datadf, sweepdf)
    return datadf.dropna()


def fill_datadf(datadf: pd.DataFrame, sweepdf: pd.DataFrame):
    sweepdf = sweepdf.copy()
    sweepdf["chem"] = np.round(sweepdf["chem"]).astype(int)
    sweepdf["food"] = np.round(sweepdf["foodp"] * 204).astype(int)
    sweepdf["dividing"] = sweepdf["tau"] - 1
    sweepdf = sweepdf.set_index(["sim", "chem", "food"])

    missing = datadf.index[datadf["dividing"].isna()]
    rounded = missing.set_levels([
        [round(i.right) for i in missing.levels[1]],
        [round(i.right) for i in missing.levels[2]]
    ], level=[1, 2])
    datadf.loc[missing, ["dividing", "bitstring"]] = sweepdf.loc[rounded, ["dividing", "bitstring"]].values


def run_pca(pcadf: pd.DataFrame):
    pca = PCA(n_components=2)
    pcomponents = pca.fit_transform(pcadf)
    resdf = pd.DataFrame(pcomponents, index=pcadf.index, columns=["pc1", "pc2"])
    resdf.index.name = "sim"
    return resdf, pca


def make_pcadf(datadf: pd.DataFrame, Jmed, Jalpha, Jweights):
    datadf = datadf.reset_index(["chemcat", "foodcat"])
    pca_dividing = datadf.pivot(columns=["chemcat", "foodcat"], values="dividing")
    gammadf = pd.merge(datadf, datadf, on="sim")
    jcc_table = contact_energy_table(Jweights)

    gammadf["gamma"] = [calculate_gamma(
        Jmed,
        Jalpha,
        cell_contact_energy(*unpack_bitstring(bit1), *unpack_bitstring(bit2), jcc_table)
    ) for bit1, bit2 in zip(gammadf["bitstring_x"], gammadf["bitstring_y"])]
    pca_gamma = gammadf.pivot(columns=["chemcat_x", "foodcat_x", "chemcat_y", "foodcat_y"], values="gamma")
    pcadf = pd.concat([pca_gamma, pca_dividing], axis=1)
    return preprocess_data(pcadf)


def unpack_bitstring(bitstring):
    js = bitstring.split("-")
    return int(js[0]), int(js[1])


def make_datadf(simdf: pd.DataFrame, chembins, foodbins):
    simdf = simdf.drop(simdf[simdf["food"] < 0].index)
    simdf["chemcat"] = pd.cut(simdf["grad_conc"], chembins)
    simdf["foodcat"] = pd.cut(simdf["food"], foodbins)
    return aggregate_data(simdf)


def concatenate_sims(celldfs: dict[str, pd.DataFrame]):
    dfs = []
    for k, celldf in celldfs.items():
        df = celldf.reset_index(drop=True)
        df["sim"] = k
        dfs.append(df)
    return pd.concat(dfs).sort_values("sim")


def make_parspacedf(datadf: pd.DataFrame):
    return pd.pivot_table(
        datadf.reset_index(),
        index="chemcat",
        columns="foodcat",
        values="sim",
        aggfunc="nunique"
    )


def preprocess_data(pcadf: pd.DataFrame):
    scaledmatrix = StandardScaler().fit_transform(pcadf)
    return pd.DataFrame(scaledmatrix, index=pcadf.index, columns=pcadf.columns).dropna(axis=1)


def aggregate_data(simdf: pd.DataFrame):
    simdf["dividing"] = simdf["tau"] == 2
    simdf["bitstring"] = simdf["jkey_dec"].astype(str) + "-" + simdf["jlock_dec"].astype(str)
    datadf = simdf.groupby(["sim", "chemcat", "foodcat"]).agg({
        "dividing": pd.Series.mean,
        "bitstring": lambda s: pd.Series.mode(s)[0]
    })
    return datadf


# This might improve performance when applied to the df before aggregate_data
def filter_matches(simdf: pd.DataFrame):
    ncounts = simdf.groupby(["chemcat", "foodcat"])["sim"].transform(pd.Series.nunique)
    nsims = simdf["sim"].max() + 1
    return simdf[ncounts == nsims].copy()


if __name__ == "__main__":
    main()
