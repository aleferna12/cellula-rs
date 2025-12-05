import logging
import argparse
import warnings
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from typing import List, Set
from random import seed, shuffle
from plotly.subplots import make_subplots
from colorir import *
from scripts.calculate_adh import *
from scripts.fileio import *

CellCluster = Set[int]
logger = logging.getLogger(__name__)
colors = px.colors.qualitative.Plotly
colors[0], colors[2] = colors[2], colors[0]


def get_parser():
    def run(args):
        if args.outputfile is None and args.adhesionfile is None:
            parser.error("either -o or -a must be set for output")

        if Path(args.datadir).is_file():
            celldf = parse_cell_data(args.datadir)
            if args.n is not None:
                celldf = reduce_data(celldf, args.n)
        else:
            t_filter = build_time_filter(get_time_points(args.datadir),
                                         start=args.start_time,
                                         end=args.end_time,
                                         n=args.n)
            celldf = parse_cell_data(args.datadir, time_filter=t_filter)
        celldf = celldf.drop(columns=gene_attrs)
        if args.outputfile is not None:
            plot_timeline(celldf, args.outputfile)
        if args.adhesionfile is not None:
            plot_adhesion(celldf,
                          args.adhesionfile,
                          jweights=np.fromstring(args.Jweights, sep=",", dtype=float),
                          top_n=args.top_n,
                          min_cluster=args.min_cluster,
                          n_pairs=args.n_pairs)

    parser = argparse.ArgumentParser(
        description="Plot information about the evolution of the simulation over time."
    )
    parser.add_argument("datadir", help="Directory containing the cell CSV files")
    parser.add_argument("-t",
                        "--start-time",
                        help="First time-step to plot (default: %(default)s)",
                        default=0,
                        type=int)
    parser.add_argument("-e",
                        "--end-time",
                        help="Last time-step to plot (default: %(default)s)",
                        default=-1,
                        type=int)
    parser.add_argument("-n",
                        help="Number of time-steps to plot (can be used to speed up plotting)",
                        type=int)
    parser.add_argument("--Jweights",
                        help="Comma-separated list of weights used to calculate gamma between bitstrings "
                             "(default: %(default)s)",
                        default="1,2,3,4,5,6,7,8")
    output = parser.add_argument_group("output arguments",
                                       "At least one of the following must be set.")
    output.add_argument("-o",
                        "--outputfile",
                        help="Output HTML or SVG file", )
    output.add_argument("-a",
                        "--adhesionfile",
                        help="Output HTML or SVG file for adhesion plots (these are slow)")
    adhesion = parser.add_argument_group("adhesion plot arguments")
    adhesion.add_argument("--top-n",
                          help="How many of the clusters with best fitness to plot (default: %(default)s)",
                          default=5,
                          type=int)
    adhesion.add_argument("--min-cluster",
                          help="Minimum size of cluster to be considered for the plot (default: %(default)s)",
                          default=50,
                          type=int)
    adhesion.add_argument("-p",
                          "--n-pairs",
                          help="Number of random cell pairs to query for adhesion per time-step "
                               "(default: %(default)s)",
                          default=50,
                          type=int)
    parser.set_defaults(run=run)
    return parser


def plot_timeline(celldf: pd.DataFrame, outputfile):
    logger.info("Making timeline plots")

    fig = make_subplots(
        3,
        2,
        shared_xaxes=True,
        subplot_titles=["Population size",
                        "Median self-gamma of cells",
                        "Total food accumulated",
                        "",
                        "Mean food per cell"]
    )
    # Plot food stuff
    mean_foods = celldf.groupby("time")["food"].mean()
    fig.add_trace(go.Scatter(
        x=mean_foods.index,
        y=mean_foods,
        line_color=colors[0],
        name="all",
        legendgroup=0,
    ), row=3, col=1)
    # Transform time into a category to make sure we get Nan values for mean if no cell with
    # a particular tau was alive at that point
    tau_foods = celldf.groupby(["tau", celldf["time"].astype("category")])["food"].agg([
        "size",
        "sum",
        "mean"
    ])
    for i, col in enumerate(tau_foods.columns, 1):
        for tau in [1, 2]:
            fig.add_trace(go.Scatter(
                x=tau_foods.loc[tau].index,
                y=tau_foods.loc[tau, col],
                line_color=colors[tau],
                name="mig" if tau == 1 else "div",
                stackgroup=None if col == "mean" else col,
                legendgroup=tau,
                showlegend=True if col == "mean" else False
            ), row=i, col=1)

    # Plot self-gammas
    med_gammas = celldf.groupby("time")["self_gamma"].median()
    fig.add_trace(go.Scatter(
        x=med_gammas.index,
        y=med_gammas,
        name="all",
        legendgroup=0,
        line_color=colors[0],
        showlegend=False
    ), row=1, col=2)
    tau_gammas = celldf.groupby(["tau", celldf["time"].astype("category")])["self_gamma"].median()
    for tau in [1, 2]:
        fig.add_trace(go.Scatter(
            x=tau_gammas.loc[tau].index,
            y=tau_gammas[tau],
            name="mig" if tau == 1 else "div",
            legendgroup=tau,
            line_color=colors[tau],
            showlegend=False
        ), row=1, col=2)
    fig.update_layout(xaxis_showticklabels=True, template="plotly_white")

    write_plot(fig, outputfile)


def plot_adhesion(celldf: pd.DataFrame,
                  outputfile,
                  jweights,
                  top_n=5,
                  min_cluster=50,
                  n_pairs=50):
    logger.info("Making adhesion plots")

    fig = make_subplots(
        3, 2,
        shared_xaxes=True,
        subplot_titles=["% of adhering cells",
                        "Median gamma of random pairs",
                        f"Top {top_n} clusters by total food (larger than {min_cluster} cells)",
                        "Median gamma of neighbouring cells"]
    )
    x = []
    multicel_perc = []
    gamma_plots = {
        "all x all": [],
        "mig x mig": [],
        "div x div": [],
        "mig x div": []
    }
    cluster_plot = {
        "time": [],
        "mean_food": [],
        "total_food": [],
        "size": []
    }
    pairs = []

    time_steps = celldf["time"].unique()
    pbar = Counter(total=len(time_steps), desc='Time-steps')
    for time_step in time_steps:
        x.append(int(time_step))
        sdf = celldf[celldf["time"] == time_step]

        clusters = get_adhering_clusters(sdf)
        unicells = 0
        total_pop = 0
        for c in clusters:
            # Notice we dont use min_cluster to access this parameter, that only determines
            # which clusters to plot as a bubble!
            if len(c) == 1:
                unicells += 1
            else:
                cluster_plot["time"].append(int(time_step)),
                cluster_food = sdf.loc[sdf["sigma"].isin(c), "food"].agg(["mean", "sum"])
                cluster_plot["mean_food"].append(cluster_food["mean"])
                cluster_plot["total_food"].append(cluster_food["sum"])
                cluster_plot["size"].append(len(c))
            total_pop += len(c)
        multicel_perc.append(1 - unicells / total_pop)

        neigh_gammas = []
        mig_mig_gammas = []
        div_div_gammas = []
        mig_div_gammas = []
        adh_table = contact_energy_table(jweights)
        it = zip(
            sdf["sigma"],
            sdf["tau"],
            sdf["Jmed"],
            sdf["neighbour_list"],
            sdf["Jneighbour_list"],
            sdf["jkey_dec"],
            sdf["jlock_dec"]
        )
        jmed, jalpha = sdf.iloc[0][["Jmed", "Jalpha"]]
        for sigma, tau, sJmed, neighs, Jneighs, key, lock in it:
            for i in np.random.randint(0, len(sdf), min(n_pairs, len(sdf))):
                cs = sdf.iloc[i]
                jcc = cell_contact_energy(key, lock, cs["jkey_dec"], cs["jlock_dec"], adh_table)
                pairs.append({
                    "time": time_step,
                    "tau1": tau,
                    "tau2": cs["tau"],
                    "gamma": calculate_gamma(jmed, jalpha, jcc)
                })

            for neigh, neighJ in zip(np.fromstring(neighs, sep=' ', dtype=int),
                                     np.fromstring(Jneighs, sep=' ')):
                if neigh == 0:
                    continue
                neigh_gamma = sJmed - neighJ / 2
                if tau == sdf.loc[(time_step, neigh), "tau"]:
                    if tau == 1:
                        mig_mig_gammas.append(neigh_gamma)
                    else:
                        div_div_gammas.append(neigh_gamma)
                else:
                    mig_div_gammas.append(neigh_gamma)
                neigh_gammas.append(neigh_gamma)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gamma_plots["all x all"].append(np.median(neigh_gammas))
            gamma_plots["mig x mig"].append(np.median(mig_mig_gammas))
            gamma_plots["div x div"].append(np.median(div_div_gammas))
            gamma_plots["mig x div"].append(np.median(mig_div_gammas))

        pbar.update()

    # Plot prevalence of multicellularity
    fig.add_trace(go.Scatter(
        x=x,
        y=multicel_perc,
        name="",
        legendgroup=0,
        showlegend=False,
    ), row=1, col=1)

    # Plot clusters
    cdf = pd.DataFrame.from_dict(cluster_plot)
    # Select cluster larger than min_cluster
    cdf = cdf[cdf["size"] >= min_cluster]
    # Select top_n values from each time-step
    topcdf = cdf.sort_values("total_food", ascending=False).groupby("time").head(top_n)
    fig.add_trace(go.Scatter(
        x=topcdf["time"],
        y=topcdf["total_food"],
        mode="markers",
        name="",
        marker=dict(color=topcdf["mean_food"],
                    colorscale="Magma_r",
                    colorbar=dict(title="Mean food per cell",
                                  len=0.3,
                                  lenmode="fraction",
                                  yanchor="bottom",
                                  y=0),
                    size=15 * topcdf["size"] / min_cluster,
                    sizemode="area"),
        customdata=np.stack((topcdf["mean_food"], topcdf["size"]), axis=1),
        hovertemplate="time: %{x}<br>total: %{y:d}<br>mean: %{customdata[0]:.2f}<br>"
                      "size: %{customdata[1]}",
        showlegend=False
    ), row=2, col=1)

    pairsdf = pd.DataFrame(pairs)
    pairsdf["time"] = pd.Categorical(pairsdf['time'], categories=x)

    def agg_pairs(df):
        return df.groupby("time")["gamma"].median()

    fig.add_trace(go.Scatter(
        x=x,
        y=agg_pairs(pairsdf),
        name="all x all",
        legendgroup=0,
        line_color=colors[0]
    ), row=1, col=2)
    fig.add_trace(go.Scatter(
        x=x,
        y=agg_pairs(pairsdf[pairsdf["tau1"] + pairsdf["tau2"] == 2]),
        name="mig x mig",
        legendgroup=1,
        line_color=colors[1]
    ), row=1, col=2)
    fig.add_trace(go.Scatter(
        x=x,
        y=agg_pairs(pairsdf[pairsdf["tau1"] + pairsdf["tau2"] == 4]),
        name="div x div",
        legendgroup=2,
        line_color=colors[2]
    ), row=1, col=2)
    fig.add_trace(go.Scatter(
        x=x,
        y=agg_pairs(pairsdf[pairsdf["tau1"] + pairsdf["tau2"] == 3]),
        name="mig x div",
        legendgroup=3,
        line_color=colors[3]
    ), row=1, col=2)

    # Plot gammas with neighbours
    for i, (name, ys) in enumerate(gamma_plots.items()):
        fig.add_trace(go.Scatter(
            x=x,
            y=ys,
            name=name,
            legendgroup=i,
            line_color=colors[i],
            showlegend=False
        ), row=2, col=2)

    fig.update_yaxes(range=[0, 1], row=1, col=1)
    fig.update_yaxes(title="total food", row=2, col=1)
    fig.update_xaxes(showticklabels=True)
    fig.update_layout(template="plotly_white")
    write_plot(fig, outputfile)


def get_adhering_neighbours(celldf: pd.DataFrame):
    """Return a list of lists with the neighbours that each cell is adhering to."""
    if not celldf["time"].nunique() == 1:
        raise ValueError("make sure this function is called for a single time point")

    adh_neigh_lists = []
    for neighs, Jneighs, Jmed in zip(celldf["neighbour_list"],
                                     celldf["Jneighbour_list"],
                                     celldf["Jmed"]):
        ad_neigh_list = []
        for neigh, neighJ in zip(np.fromstring(neighs, sep=' ', dtype=int),
                                 np.fromstring(Jneighs, sep=' ')):
            if neigh != 0 and Jmed - neighJ / 2 > 0:
                ad_neigh_list.append(neigh)
        adh_neigh_lists.append(ad_neigh_list)
    return adh_neigh_lists


def get_adhering_clusters(celldf: pd.DataFrame) -> List[CellCluster]:
    """Returns a list of list representing the different clusters of adhering cells in the
    system."""

    adh_neigh_lists = get_adhering_neighbours(celldf)
    clusters = []
    neigh_set = set(celldf["sigma"])
    neigh_dict = dict(zip(celldf["sigma"], adh_neigh_lists))
    while neigh_set:
        cluster = set()
        qneighs = {neigh_set.pop()}
        while qneighs:
            neigh = qneighs.pop()
            # Prevent infinite recursion
            if neigh not in cluster:
                cluster.add(neigh)
                qneighs.update(neigh_dict[neigh])

        neigh_set -= cluster
        clusters.append(set(cluster))

    return clusters


def get_cluster_colors(clusters: list[CellCluster], min_cluster=2, outcast_color="ffffff"):
    pal = StackPalette.load("carnival")
    grad = PolarGrad(pal, color_sys=HCLuv)
    # Only allocate colors for clusters with at least min_cluster neighbours
    n_colors = len([c for c in clusters if len(c) >= min_cluster])
    color_list = []
    # There is a bug that prevents colorir.Grad.n_colors(1)
    if n_colors == 1:
        color_list = [grad.n_colors(3)[1]]
    elif n_colors > 1:
        color_list = grad.n_colors(n_colors)

    # Want consistent results for n clusters
    seed(0)
    cluster_colors = []
    for cluster in clusters:
        if len(cluster) >= min_cluster:
            shuffle(color_list)
            color = color_list.pop()
        else:
            color = config.DEFAULT_COLOR_FORMAT.format(outcast_color)
        cluster_colors.append(color)

    return cluster_colors
