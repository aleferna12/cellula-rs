import seaborn as sns
import colorir as cl
import graphviz as gv
import matplotlib.pyplot as plt
import scipy.stats as ss
from matplotlib.collections import LineCollection, RegularPolyCollection

from scripts.data_processing import filter_kde, make_density_matrix
from scripts.calculate_adh import *
from scripts.fileio import *
from scripts.config import CS

SQRT2 = 2 ** (1 / 2)
logger = logging.getLogger(__name__)
div_grad = cl.Grad(["#ffaeae", "#cacaca", "#b0e5ff"], color_format=cl.MATPLOTLIB_COLOR_FORMAT)


def get_parser():
    def run(args):
        if Path(args.datadir).is_file():
            celldf = parse_cell_data(args.datadir)
            if args.n is not None:
                celldf = reduce_data(celldf, args.n)
        else:
            t_filter = build_time_filter(get_time_points(args.datadir),
                                         start=args.start_time,
                                         n=args.n)
            celldf = parse_cell_data(args.datadir, time_filter=t_filter, n_processes=args.n_processes)
        celldf = celldf.drop(columns=gene_attrs)
        plot_regulation(celldf,
                        plotfile=args.plotfile,
                        graphfile=args.graphfile,
                        jweights=np.fromstring(args.Jweights, sep=",", dtype=float),
                        latticesize=args.latticesize,
                        gradscale=args.gradscale,
                        gamma_thresh=args.gamma_thresh,
                        n_bins=args.n_bins,
                        plot_kernel=args.plot_kernel,
                        incl_thres=args.include_thres,
                        grid_size=args.grid_size,
                        line_width=args.line_width,
                        rel_taus=args.rel_taus,
                        discrete_taus=args.discrete_taus)
        logger.info("Finished")

    parser = argparse.ArgumentParser(
        description="Plot information about how cells regulate behaviour according to parameters."
    )
    parser.add_argument("datadir",
                        help="Directory containing the cell CSV files")
    parser.add_argument("plotfile",
                        help="Output SVG file of the main plot")
    parser.add_argument("graphfile",
                        help="Output SVG file of the network of bitstring interactions")
    parser.add_argument("latticesize",
                        help="Size of the lattice in the simulation",
                        type=int)
    parser.add_argument("-n",
                        help="Number of time-steps to plot (can be used to speed up plotting)",
                        default=None,
                        type=int)
    parser.add_argument("-t",
                        "--start-time",
                        help="First time-step to plot (default: %(default)s)",
                        default=1e6,
                        type=int
                        )
    parser.add_argument("-k",
                        "--plot-kernel",
                        help="Use to plot the kernel lines that represent cell population",
                        action="store_true"
                        )
    parser.add_argument(
        "--rel-taus",
        help="Use to plot tau as how much more likely it is that cells are dividing "
             "compared to the population mean",
        action="store_true"
    )
    parser.add_argument(
        "--discrete-taus",
        help="Use to discretize the range of taus",
        default=None,
        type=int
    )
    parser.add_argument(
        "-l",
        "--line-width",
        help="Width of the lines around isoregions",
        default=4,
        type=float
    )
    parser.add_argument(
        "-g",
        "--gradscale",
        help="Gradscale parameter used to run simulation (default: %(default)s)",
        default=5,
        type=float
    )
    parser.add_argument(
        "-b",
        "--n-bins",
        help="In how many bins the data will be divided for plotting "
             "(default: use data to decide)",
        default=None,
        type=int
    )
    parser.add_argument(
        "--include-thres",
        help="How much of the real data to include in the plot (determines the ranges of the plot), expressed "
             "in the range [0, 1], (default: %(default)s)",
        default=0.9,
        type=float
    )
    parser.add_argument(
        "--grid-size",
        help="Size of the grid used for KDE-based filtering (default: %(default)s)",
        default=100,
        type=int
    )
    parser.add_argument(
        "--gamma-thresh",
        help="Maximum gamma for an interaction to be considered weak (plotted as a dashed line) "
             "(default: %(default)s)",
        default=6.5,
        type=float
    )
    parser.add_argument(
        "--Jweights",
        help="Comma-separated list of weights used to calculate gamma between bitstrings "
             "(default: %(default)s)",
        default="1,2,3,4,5,6,7,8"
    )
    parser.add_argument("--n-processes",
                        help="How many processes are used to run the analysis"
                             "(default: %(default)s)",
                        default=1,
                        type=int)
    parser.set_defaults(run=run)
    return parser


def plot_regulation(celldf: pd.DataFrame,
                    plotfile,
                    graphfile,
                    jweights,
                    latticesize,
                    gradscale,
                    gamma_thresh,
                    n_bins=None,
                    plot_kernel=False,
                    incl_thres=1,
                    grid_size=100,
                    line_width=2,
                    discrete_taus=None,
                    rel_taus=False):
    logger.info("Plotting regulation in parameter space")
    plt.rc('font', size=14)
    celldf = celldf.sort_values("tau")
    celldf["bitstring"] = [f"{jkey}-{jlock}" for jkey, jlock in zip(celldf["jkey_dec"], celldf["jlock_dec"])]
    celldf["chem"] = [norm_chem(grad, latticesize, gradscale) for grad in celldf["grad_conc"]]
    celldf["tau_str"] = celldf["tau"].astype(str)

    cont_table = pd.crosstab(celldf["bitstring"], celldf["tau_str"])
    if cont_table.shape[0] >= 2:
        ass = ss.contingency.association(cont_table, correction=True, method='cramer')
        chitest = ss.chi2_contingency(cont_table, correction=True)
        logger.info("Association between bitstring and tau expression is: "
                    f"{ass:.4f} (p-value {chitest.pvalue:.4f})")

    # Cuts down celldf for plotting
    if incl_thres < 1:
        dm = make_density_matrix(celldf, grid_size=grid_size)
        celldf = filter_kde(celldf, include_thres=incl_thres, density_matrix=dm)

    fig, ax = plt.subplots()
    if plot_kernel:  # Mostly for debugging
        # bitstring_colors = cl.StackPalette.load("set1")[2:]
        # sns.set_palette(bitstring_colors)
        sns.set_palette(["#c72114", "#045a8d"])  # Use to plot tau_str in kde
        sns.kdeplot(
            celldf,
            x="chem",
            y="food",
            hue="tau_str",
            levels=3,
            thresh=0.01,
            legend=False,
            zorder=1,
            gridsize=grid_size,
            ax=ax
        )

    n_chem = celldf["chem"].unique().size
    n_bins = n_chem if n_bins is None else min(n_bins, n_chem)
    # Separates the bins in categories so there is not much overlap in the scatter plot
    celldf["chembin"] = pd.cut(celldf["chem"], bins=n_bins)
    celldf["foodbin"] = pd.cut(celldf["food"], bins=n_bins)
    # Plot a marker for each food and chem observed, colored by the percentage of cells dividing
    bitmode = celldf.groupby(["chembin", "foodbin"], sort=False, observed=False).agg(
        {"bitstring": lambda s: pd.Series.mode(s)[0],
         "tau": lambda s: s.mean() - 1}
    ).reset_index()
    bitstrings = bitmode["bitstring"].dropna().unique()

    heatmap_data = pd.pivot_table(
        bitmode,
        index="foodbin",
        columns="chembin",
        values="tau",
        dropna=False,
        observed=False
    )
    # Sometimes when n_bins is only slightly smaller than n_chem we end up with empty databins on the xaxis
    # This fixes the problem
    # TODO: currently it always interpolates 1 block as long as the NaN is not at the margins
    #       implement the workaround described here: https://stackoverflow.com/questions/67128364/how-to-limit-pandas-interpolation-when-there-is-more-nan-than-the-limit
    # heatmap_data = heatmap_data.interpolate(axis=1, limit=1, limit_area="inside")
    sns.heatmap(
        heatmap_data,
        cmap=div_grad.to_cmap(),
        cbar=False,
        square=True,
        xticklabels=[],
        yticklabels=[]
    )

    ax.invert_yaxis()
    ax.set_xlabel("chem. signal")
    ax.set_ylabel("metabolic reserves")
    ticks = np.linspace(0, len(heatmap_data) - 1, 5, dtype=int)
    xticklabels = [f"{heatmap_data.columns[x].mid:.1f}" for x in ticks]
    yticklabels = [f"{heatmap_data.index[y].mid:.0f}" for y in ticks]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(xticklabels)
    ax.set_yticklabels(yticklabels)
    ax.spines[:].set_visible(True)

    color_list = cl.Grad(cl.StackPalette.load("carnival")).n_colors(len(bitstrings), include_ends=True)
    bitcolors = dict(zip(bitstrings, color_list))
    iso_data = pd.pivot_table(
        bitmode,
        index="foodbin",
        columns="chembin",
        values="bitstring",
        dropna=False,
        aggfunc="first",
        observed=False
    )
    for bitstring, color in bitcolors.items():
        # Havent quite figured out where the 1.45 comes from, for now Ive just set on it empirically
        lines, corners = iso_lines(iso_data.values == bitstring, line_width / ax.bbox.width / 1.45 * n_bins)

        ax.add_collection(LineCollection(lines, colors=color, lw=line_width))
        ax.add_collection(RegularPolyCollection(
            numsides=4,
            rotation=np.radians(45),
            sizes=((line_width * 1.25) ** 2,),
            facecolors=(color,),
            offsets=corners,
            offset_transform=ax.transData,
        ))

    jmed, jalpha = celldf.iloc[0][["Jmed", "Jalpha"]]
    bit_taus = celldf.groupby("bitstring")["tau"].mean() - 1
    div_grad_loc = div_grad
    if rel_taus:
        epsilon = 1e-8
        div_mean = celldf["tau"].mean() - 1 + epsilon
        bit_taus = bit_taus.clip(lower=epsilon) / div_mean
        # How many times over the base rate
        bit_taus[:] = np.log2(bit_taus)
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

    logger.info("Bitstrings and respective % of dividing cells are:\n"
                f"{bit_taus.loc[bitstrings].reset_index().values}")

    nodecolors = {bit: div_grad_loc(x) for bit, x in bit_taus.items()}
    g = make_graph(bitstrings,
                   jmed,
                   jalpha,
                   jweights,
                   gamma_thresh,
                   colors=nodecolors,
                   bordercolors=bitcolors)
    logger.info("Writing output files")
    fig.savefig(plotfile, bbox_inches="tight")
    g.render(outfile=graphfile, cleanup=True)


def iso_lines(matrix, adjust):
    matrix = np.pad(matrix, (0, 1))
    v = np.diff(matrix.astype(int), axis=1, prepend=0).T
    h = np.diff(matrix.astype(int), axis=0, prepend=0).T

    lines = []
    for direction in [1, -1]:
        for axis in range(2):
            mat = v if axis == 0 else h
            ind_start = np.argwhere(mat == direction).astype(float)
            ind_start[:, axis] += adjust * direction
            ind_end = ind_start.copy()
            ind_end[:, 1 - axis] += 1
            lines.append(np.stack([ind_start, ind_end], axis=1))

    corners = []
    padded_v = np.pad(v[1:], ((0, 1), (0, 0)))
    padded_h = np.pad(h[:, 1:], ((0, 0), (0, 1)))
    corners.append(np.argwhere((v == -1) & (h == -1)).astype(float) + [-adjust, -adjust])
    corners.append(np.argwhere((padded_v == 1) & (h == -1)).astype(float) + [1 + adjust, -adjust])
    corners.append(np.argwhere((v == -1) & (padded_h == 1)).astype(float) + [-adjust, 1 + adjust])
    corners.append(np.argwhere((padded_v == 1) & (padded_h == 1)).astype(float) + [1 + adjust, 1 + adjust])
    return np.concatenate(lines), np.concatenate(corners)


def make_graph(bitstrings,
               jmed,
               jalpha,
               jweights,
               gamma_thresh,
               markers=None,
               colors=None,
               labels=None,
               bordercolors=None,
               sizes=None,
               fontsize="",
               name="",
               graph_class=gv.Graph):
    if markers is None:
        markers = {}
    if colors is None:
        colors = {}
    if labels is None:
        labels = {}
    if bordercolors is None:
        bordercolors = {}
    if sizes is None:
        sizes = {}

    adh_table = contact_energy_table(jweights)
    max_gamma = jmed - jalpha / 2
    adh_grad = cl.Grad(
        ["ffffff", "000000"],
        domain=[0, max_gamma],
        color_format=cl.MATPLOTLIB_COLOR_FORMAT,
        # color_coords=[gamma_thresh, max_gamma]  # Include to only show edges with gamma >
        # gamma_thresh
        # color_coords=[0, 0]  # Include to make all lines black
        # color_coords=[gamma_thresh, gamma_thresh]  # Include to only show edges with high
        # gamma and make all black
    )
    g = graph_class(name)
    g.attr(nodesep="0.25", ranksep="0.5", rankdir="TB")
    g.attr("node",
           style="filled",
           width="0.5",
           height="0.5",
           penwidth="4",
           fixedsize="true")
    g.attr("edge", fontsize="10", penwidth="2", arrowhead="none", arrowtail="none")
    for i, bitstring1 in enumerate(bitstrings):
        name1 = name + bitstring1
        color = colors.get(bitstring1, CS.black)
        size = str(sizes.get(bitstring1, 0.5))
        g.node(
            name1,
            fillcolor=color,
            fontcolor=get_fontcolor(color),
            fontsize=str(fontsize),
            shape=markers.get(bitstring1, "circle"),
            label=labels.get(bitstring1, ""),
            color=bordercolors.get(bitstring1, "#ffffff"),
            width=size,
            height=size
        )
        k1, l1 = bitstring1.strip("_").split("-")
        for bitstring2 in bitstrings[i:]:
            name2 = name + bitstring2
            k2, l2 = bitstring2.strip("_").split("-")
            jcc = cell_contact_energy(int(k1), int(l1), int(k2), int(l2), adh_table)
            edge_gamma = calculate_gamma(jmed, jalpha, jcc)
            if edge_gamma < 0:
                continue

            style = "solid" if edge_gamma > gamma_thresh else "dashed"
            g.edge(name1,
                   name2,
                   color=adh_grad(edge_gamma),
                   # label=str(edge_gamma),  # Uncomment to show gamma as text
                   style=style)
    return g


def get_fontcolor(color):
    if cl.config.DEFAULT_COLOR_FORMAT.format(color).cielab().l < 50:
        return CS.white
    else:
        return CS.black


def chem_to_dist(chem, lattsize, gradscale):
    lattdiag = SQRT2 * lattsize
    return lattdiag - 100 * (chem - 1) / gradscale


def norm_chem(chem, lattsize, gradscale, norm=100):
    lattdiag = SQRT2 * lattsize
    return norm * (chem - 1) / (gradscale * lattdiag / 100)
