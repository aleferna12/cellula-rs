import numbers
from tempfile import TemporaryFile, NamedTemporaryFile, TemporaryDirectory

import numpy as np
import pandas as pd
import networkx as nx
import plotly.graph_objects as go
import colorir as cl
import graphviz as gv
from functools import reduce
from enlighten import Counter
from statsmodels.stats.weightstats import DescrStatsW

from scripts.analyses.make_netgraphs import count_mut
from scripts.data_processing import make_density_matrix, filter_kde, get_parameter_range
from scripts.data_processing import logger as data_logger
from scripts.fileio import *
from scripts.analyses.make_trees import make_trees, get_longest_trees
from scripts.sweep import parse_sweep_args, sweep_cell, add_inestimable_sweep_args
from scripts.calculate_adh import *
from scripts.analyses.plot_regulation import make_graph, div_grad
from scripts.config import CS

logger = logging.getLogger(__name__)


def get_parser():
    def run(args):
        data_logger.setLevel(logging.WARNING)

        if args.ancfile is None and args.bitfile is None:
            parser.error("either '--ancfile' or '--bitfile' should be set for output")

        celldf = parse_cell_data(args.datadir, n_processes=args.n_processes)
        celldf = celldf[celldf["time"] != 0]
        sweepkwargs = parse_sweep_args(args)
        jmed, jalpha = celldf.iloc[0][["Jmed", "Jalpha"]]
        jweights = np.fromstring(args.Jweights, sep=",", dtype=float)
        ancdf = get_ancestors(celldf, args.n_times)
        anc_interval = ancdf.iloc[1]["time"] - ancdf.iloc[0]["time"]
        logger.info(f"Ancestors are: {anc_interval} MCS apart")

        anc_bitfreqs = {}
        pop_bitfreqs = {}
        anc_sweeps = {}
        pbar = enlighten.Counter(desc="Ancestral populations analysed", total=len(ancdf))
        for (time, _), row in ancdf.iterrows():
            popdf = celldf[celldf["time"].between(
                time - args.time_span,
                time,
                inclusive="both"
            )]
            dm = make_density_matrix(popdf, grid_size=100, n_processes=args.n_processes)
            pr = get_parameter_range(filter_kde(
                popdf,
                include_thres=args.include_thres,
                density_matrix=dm
            ))
            anc_sweeps[time] = sweep_cell(
                row,
                parameter_range=pr,
                **sweepkwargs
            )
            anc_bitfreqs[time] = estimate_bitfreqs(anc_sweeps[time], dm)
            pop_bitfreqs[time] = calculate_bifreqs(popdf)
            pbar.update()

        ancdatadf = make_datadf(
            anc_bitfreqs,
            anc_bitfreqs if args.self_gamma else pop_bitfreqs,
            jmed,
            jalpha,
            jweights
        )

        if args.ancfile is not None:
            plot_ancestry(ancdatadf, args.ancfile, args.div_gamma)

        if args.bitfile is not None:
            plot_bitstring_evolution(
                ancdatadf["strat"],
                anc_sweeps,
                anc_bitfreqs,
                pop_bitfreqs,
                args.bitfile,
                jmed,
                jalpha,
                jweights,
                args.gamma_thresh,
                minsize=args.minsize,
                maxsize=args.maxsize,
                pop_bits=args.pop_bits,
                numerate_bitstrings=args.numerate,
                rel_taus=args.rel_taus,
                discrete_taus=args.discrete_taus,
                dev_overlap=args.dev_overlap,
                hamm_thresh=args.hamm_thresh,
                first_ancestor=args.first_ancestor
            )
        logger.info("Finished")

    parser = argparse.ArgumentParser(
        description="Plot data about the ancestors of a simulation."
    )
    parser.add_argument("datadir",
                        help="Directory containing the cell CSV files")
    output = parser.add_argument_group("output arguments",
                                       description="What files to output. At least one must be set.")
    output.add_argument("-a",
                        "--ancfile",
                        help="Output HTML or SVG file for the evolution of the ancestral strategy")
    output.add_argument("-b",
                        "--bitfile",
                        help="Output SVG file for the evolution of bitstring networks")
    parser.add_argument("-t",
                        "--n-times",
                        help="How many time points will be sampled (default: %(default)s)",
                        default=10,
                        type=int)
    parser.add_argument("-g",
                        "--gamma-thresh",
                        help="Gamma threshold for dashing a line connecting two bitstrings "
                             "(default: %(default)s)",
                        default=6.5,
                        type=float)
    parser.add_argument("-s",
                        "--time-span",
                        help="Length of the time period where the population is considered to not have changed, "
                             "the more data points you have the smaller this can be "
                             "(default: %(default)s)",
                        default=2e5,
                        type=float)
    parser.add_argument("--minsize",
                        help="Minimum size for a node when its expressed by 0.%% of the pop. (default: %(default)s)",
                        default=0.5,
                        type=float)
    parser.add_argument("--maxsize",
                        help="Maximum size for a node when its expressed by 100.%% of the pop. (default: %(default)s)",
                        default=0.5,
                        type=float)
    parser.add_argument("--Jweights",
                        help="Comma-separated list of weights used to calculate gamma between bitstrings "
                             "(default: %(default)s)",
                        default="1,2,3,4,5,6,7,8")
    parser.add_argument("--n-processes",
                        help="How many processes are used to run the analysis"
                             "(default: %(default)s)",
                        default=1,
                        type=int)
    parser.add_argument("--pop-bits",
                        help="How many bitstrings to show from the population (per cell type)"
                             "(default: %(default)s)",
                        default=0,
                        type=int)
    parser.add_argument("--numerate",
                        help="If set, bitstring labels are numerated",
                        action="store_true")
    parser.add_argument("--div-gamma",
                        help="If set, combine migmig gamma and migdiv gamma into div gamma",
                        action="store_true")
    parser.add_argument("--self-gamma",
                        help="If set, uses the ancestor gamma with itself for the plot",
                        action="store_true")
    parser.add_argument(
        "--rel-taus",
        help="Use to plot tau as how much more likely it is that cells are dividing "
             "compared to the population mean",
        action="store_true"
    )
    parser.add_argument(
        "--dev-overlap",
        help="Use to indicate that a bitstring co-option event requires that "
             "the co-opted and evolved bitstrings overlap in developmental space",
        action="store_true"
    )
    parser.add_argument(
        "--hamm-thresh",
        help="Maximum hamming distance between bitstrings where they are"
             " still considered homologous (default: %(default)s)",
        type=int,
        default=4
    )
    parser.add_argument(
        "--first-ancestor",
        help="Use to indicate that the cooption algorithm should track all the way to the first ever ancestor "
             "(default behaviour is stopping at the last ancestor that expresses the ancestral life cycle)",
        action="store_true"
    )
    parser.add_argument(
        "--discrete-taus",
        help="Use to discretize the range of taus",
        default=None,
        type=int
    )
    add_inestimable_sweep_args(parser)
    parser.set_defaults(run=run)
    return parser


def plot_ancestry(datadf: pd.DataFrame, outfile, div_gamma):
    fig = go.Figure()

    # Draw bg recs
    strat_colors = cl.Palette.load("strats").blend(CS.white, 0.5)
    strat_grad = cl.Grad(
        strat_colors.colors[:-1],
        color_coords=[1 / 8, 3 / 8, 5 / 8, 7 / 8],
        discrete=True
    )
    times = datadf.index.get_level_values(0) - 1e8
    points = list(zip(times, datadf["strat"], datadf["div_gamma"]))
    points.append(points[-1])
    x0 = 0
    # Each iteration plots the strategy of the previous ancestor
    for i in range(1, len(points)):
        time, strat, d_gamma = points[i - 1]
        next_time, next_strat, _ = points[i]
        x1 = (time + next_time) / 2
        fig.add_vrect(
            x0=x0,
            x1=x1 if strat != next_strat else x1 * 1.01,  # Gets rid of visual glitch
            fillcolor=get_strat_color(strat, strat_colors, d_gamma),
            line_width=0,
            opacity=1,
            layer="below",
        )
        x0 = x1

    fig.add_traces([
        go.Scatter(
            x=times,
            y=datadf["migmig_gamma"],
            mode="lines+markers",
            line_color=CS.redlipstick,
            name="mig-mig"
        ),
        go.Scatter(  # Dummy trace for the color bar
            x=[None],
            y=[None],
            showlegend=False,
            marker=go.scatter.Marker(
                colorscale=strat_grad.to_plotly_colorscale(),
                showscale=True,
                cmin=0,
                cmax=1,
                colorbar=go.scatter.marker.ColorBar(
                    tickvals=strat_grad.color_coords,
                    ticktext=["undefined", "unicellular", "multicellular", "mixed"],
                    tickfont=go.scatter.marker.colorbar.Tickfont(size=14),
                    len=0.5,
                    yanchor="bottom",
                    y=0
                )
            )
        )
    ])
    if div_gamma:
        fig.add_trace(go.Scatter(
            x=times,
            y=datadf["div_gamma"],
            mode="lines+markers",
            line_color=cl.Hex("#005199"),
            name="div"
        ))
    else:
        fig.add_traces([
            go.Scatter(
                x=times,
                y=datadf["migdiv_gamma"],
                mode="lines+markers",
                line_color=CS.petunia,
                name="mig-div"
            ),
            go.Scatter(
                x=times,
                y=datadf["divdiv_gamma"],
                mode="lines+markers",
                line_color=CS.deeplagoon,
                name="div-div"
            )])
    fig.update_layout(width=600, height=300)
    fig.update_xaxes(title="time", range=[times.min(), times.max()])
    fig.update_yaxes(title="gamma", zeroline=False, range=[-19, 19])
    write_plot(fig, outfile)


def get_strat_color(strat, strat_colors, div_gamma):
    mul_grad = cl.Grad(strat_colors[["mul_propagules", "multicellular"]], domain=[0, 18])
    return strat_colors[strat] if strat != "multicellular" else mul_grad(max(0, div_gamma))


def bitstring_graph_properties(bitfreqs, minsize, maxsize, rel_taus=False):
    bit_div = {}
    bit_sizes = {}
    for i in bitfreqs:
        time_bf = bitfreqs[i]
        rel_freq = time_bf.sum(axis=1)
        rel_tau = time_bf.div(rel_freq, axis=0)["div"]
        if rel_taus:
            epsilon = 1e-8
            tot_tau = time_bf.sum(axis=0)
            rel_tau = rel_tau.clip(lower=epsilon) / (tot_tau["div"] + epsilon)
            # How many times over the base rate
            rel_tau[:] = np.log2(rel_tau)
        bit_div[i] = rel_tau.to_dict()
        # TODO: represent this in another way? Size is a bit weird
        if minsize == maxsize:
            bit_sizes[i] = dict(zip(time_bf.index, [minsize] * len(time_bf)))
        else:
            div = 1 / (maxsize - minsize)
            sum_ = minsize * div
            bit_sizes[i] = ((time_bf.div(time_bf.sum(axis=0)).max(axis=1) + sum_) / div).to_dict()
    return bit_div, bit_sizes


def plot_bitstring_evolution(
        anc_strats,
        anc_sweeps,
        anc_bitfreqs,
        pop_bitfreqs,
        outfile,
        jmed,
        jalpha,
        jweights,
        gamma_thresh,
        minsize=0.5,
        maxsize=0.5,
        pop_bits=0,
        numerate_bitstrings=False,
        rel_taus=False,
        discrete_taus=None,
        dev_overlap=False,
        hamm_thresh=4,
        first_ancestor=False
):
    for i, pop_bf in pop_bitfreqs.items():
        mig = pop_bf["mig"].sort_values(ascending=False)[:pop_bits]
        div = pop_bf["div"].sort_values(ascending=False)[:pop_bits]
        filtered = pop_bf.loc[list(set(mig.index) | set(div.index))]
        pop_bitfreqs[i] = filtered / filtered.values.sum()

    anc_div, anc_size = bitstring_graph_properties(anc_bitfreqs, minsize, maxsize, rel_taus)
    pop_div, pop_size = bitstring_graph_properties(pop_bitfreqs, minsize, maxsize, rel_taus)

    bit_div = {}
    bit_size = {}
    bit_markers = {}
    for i, div in pop_div.items():
        bit_div[i] = anc_div[i].copy()
        bit_size[i] = anc_size[i].copy()
        bit_markers[i] = {}
        size = pop_size[i]
        for bit in div:
            bit_div[i][bit + "_"] = div[bit]
            bit_size[i][bit + "_"] = size[bit]
            bit_markers[i][bit + "_"] = "square"

    bitset = reduce(lambda t1, t2: pd.unique(pd.Series(list(t1) + list(t2))), anc_div.values())
    labels = {bit: str(i) if numerate_bitstrings else bit for i, bit in enumerate(bitset)}
    nbits = len(bitset)
    popbitset = reduce(lambda t1, t2: pd.unique(pd.Series(list(t1) + list(t2))), pop_div.values())
    for bit in popbitset:
        if not numerate_bitstrings or bit in labels:
            labels[bit + "_"] = labels.get(bit, bit)
        else:
            labels[bit + "_"] = str(nbits)
            nbits += 1

    div_grad_loc = div_grad
    if rel_taus:
        div_grad_loc.domain = [-2, 2]
        div_grad_loc.color_coords = [-2, 0, 2]
    if discrete_taus is not None:
        div_grad_loc.colors = [x.hex() for x in div_grad_loc.n_colors(discrete_taus)]
        div_grad_loc.color_coords = np.linspace(
            div_grad_loc.domain[0],
            div_grad_loc.domain[1],
            discrete_taus
        )
        div_grad_loc.discrete = True

    parent = gv.Digraph()
    parent.attr(compound="true")

    gs = {}
    for i, divs in bit_div.items():
        colors = {b: div_grad_loc(d) for b, d in divs.items()}
        gs[i] = make_graph(list(colors),
                           jmed,
                           jalpha,
                           jweights,
                           gamma_thresh,
                           labels=labels,
                           colors=colors,
                           markers=bit_markers[i],
                           sizes=bit_size[i],
                           name=f"cluster_{i}",
                           graph_class=gv.Digraph,
                           fontsize=8)

    for i, g in gs.items():
        g.attr(label=f"time: {i}")
        parent.subgraph(g)

    # Creates a bunch of graphs, each tracking the evolution of one of the bitstrings of the
    # last ancestor that made unicellular propagules in the dataset all the way
    # to the first ancestor
    # We have to restrict this analysis to uniprops because its hard to determine if a bitstring is responsible
    # for making multicellular propagules (gamma_mig,div is often positive)
    uniprops = np.argwhere(anc_strats == "uni_propagules").flatten()
    if len(uniprops) > 0:
        if rel_taus:
            cooption_graphs, coopdf = make_cooption_graphs(
                anc_strats.index[uniprops[-1]],
                anc_bitfreqs,
                anc_strats,
                anc_sweeps,
                anc_div,
                jmed,
                jalpha,
                jweights,
                dev_overlap,
                hamm_thresh,
                first_ancestor,
            )
            for bit, g in cooption_graphs.items():
                g.attr(label=f"bitstring: {bit}")
                parent.subgraph(g)

            coopdf.to_csv(outfile.replace(".svg", ".csv").replace(".html", ".csv"))
        else:
            logger.warning("Not including co-option graphs since rel_tau is 'False'")

    logger.info(f"Writing plot to {outfile}")
    parent.render(outfile=outfile, cleanup=True)


def make_cooption_graphs(
        last_time,
        anc_bitfreqs,
        anc_strats,
        anc_sweeps,
        anc_div,
        jmed,
        jalpha,
        jweights,
        dev_overlap,
        hamm_thresh,
        first_ancestor,
):
    last_bits = anc_bitfreqs[last_time].index
    times = np.sort(list(anc_bitfreqs.keys()))[::-1]
    times = times[times <= last_time]
    bitevols = {}
    bitevol_graphs = {}
    for bit in last_bits:
        current = bit
        bitevols[bit] = {}
        bitevol_graphs[bit] = gv.Digraph("cluster_" + bit)
        bitevol_graphs[bit].attr(
            "node",
            nodesep="0.25",
            width="0.5",
            height="0.5",
            color="#000000",
            fixedsize="true",
            fontsize="8"
        )
        for time, prev_time in zip(times, times[1:]):
            if current != "-" and current not in anc_bitfreqs[time].index:
                raise Exception("wrong current assigned")

            precursor_candidates = sorted(anc_bitfreqs[prev_time].index, key=lambda x: hamming_distance(x, bit))
            precursor, fail = select_precursor(
                precursor_candidates,
                current,
                anc_sweeps[prev_time],
                anc_sweeps[time],
                dev_overlap,
                hamm_thresh
            )

            bitevols[bit][prev_time] = (precursor, fail)

            curr_name = bit + current + str(time)
            precursor_name = bit + precursor + str(prev_time)
            bitevol_graphs[bit].node(curr_name, label=current)
            bitevol_graphs[bit].node(precursor_name, label=precursor)
            bitevol_graphs[bit].edge(
                precursor_name,
                curr_name,
                label=str(time)
            )
            current = precursor

    anc_strat = anc_strats[anc_strats != "undefined"].iloc[0]
    strat_changed = np.any((anc_strats != anc_strat) & (anc_strats != "undefined"))
    if not strat_changed:
        return bitevol_graphs

    last_anc_strat = anc_strats.index[0]
    for i, strat in anc_strats.items():
        if strat == "undefined":
            continue
        if first_ancestor or strat != anc_strat:
            break
        last_anc_strat = i
    adh_table = contact_energy_table(jweights)
    coopted = {}
    coop_info = []
    most_expressed_valid = None
    most_expression_valid = None
    for bit in last_bits:
        bit_is_mig = anc_div[last_time][bit] < 0
        if anc_strat == "unicellular" and bit_is_mig:
            continue
        if anc_strat == "multicellular" and not bit_is_mig:
            continue

        for time, (precursor, fail) in bitevols[bit].items():
            if time < last_anc_strat:
                break
            if precursor == "-" or not was_coopted(precursor, anc_strat, jmed, jalpha, adh_table):
                coopted[bit] = (time, precursor, False)
                coop_info.append((
                    bit,
                    precursor,
                    time,
                    fail if fail is not None else "func"
                ))
                break
            coopted[bit] = (time, precursor, True)
            coop_info.append((
                bit,
                precursor,
                time,
                None
            ))

        this_expression = anc_bitfreqs[last_time].loc[bit].sum()
        if most_expression_valid is None or this_expression > most_expression_valid:
            self_gamma = get_self_gamma(bit, jmed, jalpha, adh_table)
            if ((anc_strat == "unicellular" and self_gamma <= 0)
                    or (anc_strat == "multicellular" and self_gamma > 0)):
                most_expressed_valid = bit
                most_expression_valid = this_expression

    for bit, (time, precursor, bit_coopted) in coopted.items():
        if bit_coopted:
            if bit == most_expressed_valid:
                edge_color = "#0000ff"
            else:
                edge_color = "#00ff00"
        else:
            edge_color = "#ff0000"

        bitevol_graphs[bit].edge(
            bit + precursor + str(time),
            bit + bit + str(last_time),
            color=edge_color
        )

    coopdf = pd.DataFrame.from_records(
        coop_info,
        columns=["bit", "precursor", "time", "fail"]
    ).set_index("bit")
    coopdf.loc[:, "most_expressed"] = (coopdf.index == most_expressed_valid)
    return bitevol_graphs, coopdf


def get_self_gamma(bit, jmed, jalpha, adh_table):
    jkey, jlock = [int(j) for j in bit.split("-")]
    return calculate_gamma(jmed, jalpha, cell_contact_energy(jkey, jlock, jkey, jlock, adh_table))


# Assumes that the derived strategy is propagule forming
def was_coopted(precursor, anc_strategy, jmed, jalpha, adh_table):
    if precursor == "-":
        return False

    self_gamma = get_self_gamma(precursor, jmed, jalpha, adh_table)
    if (anc_strategy == "unicellular" and self_gamma <= 0) or (anc_strategy == "multicellular" and self_gamma > 0):
        return True
    return False


def select_precursor(candidates, target, prev_sweep, current_sweep, dev_overlap, hamm_thresh):
    if not candidates:
        return "-", "empty"
    elif target in candidates:
        return target, None
    else:
        candidates = hamm_select(target, candidates, hamm_thresh)
        if not candidates:
            return "-", "hamm"
        if not dev_overlap:
            return candidates[0], None
        else:
            prev_sweep = prev_sweep.set_index(
                pd.MultiIndex.from_product([
                    interval_from_left(prev_sweep["chem"].unique()),
                    interval_from_left(prev_sweep["foodp"].unique())
                ])
            )
            sweep = current_sweep.set_index(["chem", "foodp"])
            locs = sweep[sweep["bitstring"] == target].index
            locs = [loc for loc in locs if loc in prev_sweep.index]
            cands_in_sweep = prev_sweep.loc[locs, "bitstring"].unique()

            for candidate in candidates:
                if candidate in cands_in_sweep:
                    return candidate, None
            return "-", "dev"


def hamm_select(target, candidates, thresh):
    if target == "-":
        return []
    return [c for c in candidates if c != "-" and hamming_distance(c, target) <= thresh]


def interval_from_left(data):
    interval = data[1] - data[0]
    return pd.interval_range(data[0], data[-1] + interval, periods=len(data))


def get_mrca(celldf, stop_mrca=True):
    trees = make_trees(celldf, collapse_branches=False, stop_mrca=stop_mrca)
    root = get_longest_trees(trees)[0][0]
    return root.get_common_ancestor(root.get_leaves())


def get_mrcass(celldf):
    mrca = get_mrca(celldf)
    return celldf.loc[(mrca.time, int(mrca.name))]


def get_ancestors(celldf: pd.DataFrame, n_ancestors=None):
    mrca = get_mrca(celldf, False)
    u_times = np.unique(celldf["time"])
    times = u_times if n_ancestors is None else interspaced_elements(u_times[u_times <= mrca.time], n_ancestors)
    anc_indexes = []
    stop = False
    while not stop:
        if n_ancestors is None or mrca.time in times:
            anc_indexes.append((mrca.time, int(mrca.name)))
        stop = mrca.is_root()  # Do it once more for the root
        mrca = mrca.up
    return celldf.loc[anc_indexes[::-1]]


def make_datadf(bitfreqs1, bitfreqs2, jmed, jalpha, jweights):
    if set(bitfreqs1.keys()) != set(bitfreqs2.keys()):
        raise ValueError("indexes of dictionaries containing bifreqs must match")

    data = []
    adh_table = contact_energy_table(jweights)
    pbar = Counter(desc="Strategies classified", total=len(bitfreqs1))
    for index in bitfreqs1:
        bf1, bf2 = bitfreqs1[index], bitfreqs2[index]
        gamma_distr = gamma_distributions(bf1, bf2, jmed, jalpha, adh_table, max_0=False)
        datarow = {
            "migmig_gamma": get_migmig_gamma(gamma_distr),
            "divdiv_gamma": get_divdiv_gamma(gamma_distr),
            "migdiv_gamma": get_migdiv_gamma(gamma_distr),
            "div_gamma": get_div_gamma(gamma_distr),
            "strat": categorize_strategy(gamma_distr)
        }
        data.append(datarow)
        pbar.update()
    datadf = pd.DataFrame(data)
    datadf["strat"] = pd.Categorical(
        datadf["strat"],
        ordered=True,
        categories=[
            "multicellular",
            "uni_propagules",
            "unicellular",
            "undefined"
        ]
    )
    datadf.index = bitfreqs1.keys()
    pbar.close()
    return datadf


def calculate_bifreqs(celldf):
    celldf = celldf.copy()
    celldf["bitstring"] = celldf["jkey_dec"].astype(str) + "-" + celldf["jlock_dec"].astype(str)
    celldf["tau"] = pd.Categorical(celldf["tau"], [1, 2])
    return pd.pivot_table(
        celldf[["bitstring", "tau"]].value_counts(normalize=True).reset_index(),
        index="bitstring",
        columns="tau",
        values="proportion",
        fill_value=0,
        dropna=False,
        observed=False
    ).rename(columns={1: "mig", 2: "div"})


# We then can multiply the estimated bitfreqs with the real bitfreqs or something
def estimate_bitfreqs(sweepdf, density_matrix):
    food = sweepdf["foodp"] * 204  # Assumes gradstep, lattice size etc
    chem_i = density_matrix.index.get_level_values(0)
    food_i = density_matrix.index.get_level_values(1)
    sweepdf = sweepdf.loc[
        (sweepdf["chem"].between(chem_i[0].left, chem_i[-1].right)) &
        (food.between(food_i[0].left, food_i[-1].right))
        ].copy()

    sweepdf["tau"] = pd.Categorical(sweepdf["tau"], [1, 2])
    sweepdf["pdf"] = density_matrix.loc[zip(sweepdf["chem"], food)].values
    bitfreqs = pd.pivot_table(
        sweepdf,
        index="bitstring",
        columns="tau",
        values="pdf",
        aggfunc="sum",
        fill_value=0,
        dropna=False,
        observed=False
    ).rename(columns={1: "mig", 2: "div"})
    return bitfreqs / bitfreqs.values.sum()


def gamma_distributions(bitfreqs1, bitfreqs2, jmed, jalpha, adh_table, max_0=True):
    """Statistical approach to calculating the gamma interactions in a population.

    Bitfreqs can be either the frequencies of bitstrings in a population (calculated with 'calculate_bitfreqs') or
    the likelihood of a single cell expressing each bitstring (estimated with 'estimate_bitfreqs').
    """

    if round(bitfreqs1.values.sum() + bitfreqs2.values.sum(), 3) != 2:
        raise ValueError("Sum of frequencies of one of the data frames != 1")

    rel_gammas = {
        "gamma": [],
        "migmig": [],
        "divdiv": [],
        "migdiv": []
    }
    for tup1 in bitfreqs1.itertuples():
        k1, l1 = np.array(tup1.Index.split("-"), dtype=int)
        for tup2 in bitfreqs2.itertuples():
            k2, l2 = np.array(tup2.Index.split("-"), dtype=int)
            g = calculate_gamma(jmed, jalpha, cell_contact_energy(k1, l1, k2, l2, adh_table))
            if max_0:
                g = max(0, g)
            rel_gammas["gamma"].append(g)
            rel_gammas["migmig"].append(tup1.mig * tup2.mig)
            rel_gammas["divdiv"].append(tup1.div * tup2.div)
            rel_gammas["migdiv"].append(tup1.mig * tup2.div + tup1.div * tup2.mig)
    return pd.DataFrame(rel_gammas)


def get_divdiv_gamma(gamma_distr):
    if (gamma_distr["divdiv"] == 0).all():
        return None

    return weighted_quantiles(
        gamma_distr["gamma"],
        0.5,
        gamma_distr["divdiv"]
    )


def get_migmig_gamma(gamma_distr):
    if (gamma_distr["migmig"] == 0).all():
        return None
    return weighted_quantiles(
        gamma_distr["gamma"],
        0.5,
        gamma_distr["migmig"]
    )


def get_migdiv_gamma(gamma_distr):
    if (gamma_distr["migdiv"] == 0).all():
        return None
    return weighted_quantiles(
        gamma_distr["gamma"],
        0.5,
        gamma_distr["migdiv"]
    )


def get_div_gamma(gamma_distr):
    weights = gamma_distr["migdiv"] + gamma_distr["divdiv"]
    if (weights == 0).all():
        return None
    return weighted_quantiles(
        gamma_distr["gamma"],
        0.5,
        weights
    )

def weighted_quantiles(values, quantiles, weights):
    res = DescrStatsW(values, weights=weights).quantile(quantiles, return_pandas=False)
    return res[0] if isinstance(quantiles, numbers.Number) else res


# This was taken from the internet
def weighted_quantiles_old(values, quantiles, sample_weight=None,
                           values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    wq = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        wq -= wq[0]
        wq /= wq[-1]
    else:
        wq /= np.sum(sample_weight)
    return np.interp(quantiles, wq, values)


def categorize_strategy(gamma_distr):
    if (gamma_distr["migmig"] == 0).all() or (gamma_distr["divdiv"] == 0).all():
        return "undefined"

    migmig = get_migmig_gamma(gamma_distr)
    if migmig <= 0:
        return "unicellular"
    divdiv = get_divdiv_gamma(gamma_distr)
    if divdiv <= 0:
        return "uni_propagules"
    return "multicellular"
