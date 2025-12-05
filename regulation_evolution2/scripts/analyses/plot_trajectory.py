import matplotlib.pyplot as plt
from colorir import *
from scripts.data_processing import filter_kde
from scripts.analyses.plot_regulation import norm_chem
from scripts.fileio import *

logger = logging.getLogger(__name__)


def get_parser():
    def run(args):
        tf = build_time_filter(get_time_points(args.datadir), start=1000000)
        celldf = parse_cell_data(args.datadir, n_processes=5, time_filter=tf)
        celldf = filter_kde(celldf, 0.95)
        datadf = make_datadf(celldf, args.n_bins, args.n_bins, args.latticesize, args.gradscale)
        vecdf = make_vecdf(datadf)
        g = Grad(["#ffaeae", "#cacaca", "#b0e5ff"])
        style = ["dashed" if args.dashed_propagules and x else "solid" for x in vecdf["propagule"]]
        plt.quiver(
            vecdf["chemcat"],
            vecdf["foodcat"],
            vecdf["dchem"],
            vecdf["dfood"],
            width=args.arrow_width,
            color=[g(x) for x in vecdf["dividing"]],
            ls=style,
            edgecolor="gray",
            linewidth=1,
            angles="xy",
            scale=args.arrow_scale,  # Changes the length of arrows
            scale_units="xy"
        )
        plt.gcf().set_size_inches(5, 5)
        plt.savefig(args.outfile)
        logger.info("Finished")

    parser = argparse.ArgumentParser(
        description="Plot data about the frequency of strategies in a simulation."
    )
    parser.add_argument("datadir",
                        help="Directory containing the cell CSV files")
    parser.add_argument("outfile",
                        help="Output HTML or SVG file")
    parser.add_argument("latticesize",
                        help="Size of the lattice in the simulation",
                        type=int)
    parser.add_argument("gradscale",
                        help="Value of the 'gradscale' parameter used in the simulation",
                        type=float)
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
    parser.add_argument(
        "-b",
        "--n-bins",
        help="In how many bins the data will be divided for plotting "
             "(default: use data to decide)",
        default=None,
        type=int
    )
    parser.add_argument("--arrow-width",
                        help="Width of the vector arrows (default: %(default)s)",
                        default=0.008,
                        type=float
                        )
    parser.add_argument("--arrow-scale",
                        help="Length scale of the vector arrows (default: %(default)s)",
                        default=2,
                        type=float
                        )
    parser.add_argument("-d",
                        "--dashed-propagules",
                        help="Use to show transitions from group to unicellular as dashed lines",
                        action="store_true"
                        )
    parser.set_defaults(run=run)
    return parser


def make_datadf(celldf: pd.DataFrame, chembins, foodbins, latt, gradscale):
    s_time = np.unique(celldf["time"])
    interval = s_time[1] - s_time[0]
    logger.info(f"Predicted time interval: {interval}")
    try:
        if not np.all(s_time == np.arange(celldf["time"].min(), celldf["time"].max() + interval, interval)):
            logger.error(f"Time-series is inconsistent with predicted interval {interval}")
            logger.error(s_time)
    except ValueError:
        pass

    celldf["chem"] = [norm_chem(x, latt, gradscale) for x in celldf["grad_conc"]]
    celldf["next_time"] = celldf["time"] + interval
    celldf["unicellular"] = celldf["neighbour_list"] == "0"
    celldf = celldf[["sigma", "ancestor", "time", "next_time", "chem", "food", "tau", "unicellular"]]
    datadf = pd.merge(celldf, celldf, left_on=["ancestor", "time"], right_on=["sigma", "next_time"])

    datadf["dchem"] = datadf["chem_x"] - datadf["chem_y"]
    datadf["dfood"] = datadf["food_x"] - datadf["food_y"]
    datadf["vec_mag"] = np.sqrt(datadf["dchem"] ** 2 + datadf["dfood"] ** 2)
    datadf["uchem"] = datadf["dchem"] / datadf["vec_mag"]
    datadf["ufood"] = datadf["dfood"] / datadf["vec_mag"]
    datadf["chemcat"] = [(x.left + x.right) / 2 for x in pd.cut(
        datadf["chem_y"],
        min(chembins, datadf["chem_y"].nunique()) - 1,
        right=False
    )]
    datadf["foodcat"] = [(x.left + x.right) / 2 for x in pd.cut(
        datadf["food_y"],
        min(foodbins, datadf["chem_y"].nunique()) - 1,
        right=False
    )]
    datadf["dividing"] = datadf["tau_y"] == 2
    return datadf


def make_vecdf(datadf: pd.DataFrame):
    datadf["propagule"] = datadf["unicellular_y"] & (~datadf["unicellular_x"]) & datadf["dividing"]
    return datadf.groupby(["chemcat", "foodcat", "propagule"], observed=True, as_index=False).mean()
