import plotly.graph_objects as go
import colorir as cl

from scripts.analyses.plot_ancestry import make_datadf, estimate_bitfreqs, calculate_bifreqs, get_ancestors, \
    get_strat_color
from scripts.data_processing import make_density_matrix, filter_kde, get_parameter_range
from scripts.data_processing import logger as data_logger
from scripts.fileio import *
from scripts.sweep import sweep_cells, parse_sweep_args, add_inestimable_sweep_args, sweep_cell
from scripts.calculate_adh import *

logger = logging.getLogger(__name__)


def get_parser():
    def run(args):
        data_logger.setLevel(logging.WARNING)

        celldf = parse_cell_data(args.datadir, n_processes=args.n_processes)
        celldf = celldf[celldf["time"] != 0]
        jmed, jalpha = celldf.iloc[0][["Jmed", "Jalpha"]]
        jweights = np.fromstring(args.Jweights, sep=",", dtype=float)
        sweepkwargs = parse_sweep_args(args)

        ancdf = get_ancestors(celldf)
        if args.stop_anc:
            times = interspaced_elements(ancdf["time"].values, args.n_times)
        else:
            times = interspaced_elements(np.unique(celldf["time"].values), args.n_times)
        anc_bitfreqs = {}
        anc_pop_bitfreqs = {}
        sample_bitfreqs = {}
        pop_repeated_bitfreqs = {}
        pbar = enlighten.Counter(desc="Populations sampled", total=args.n_times)
        for time in times:
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
            pop_bitfreqs = calculate_bifreqs(popdf)

            sampledf = popdf.sample(min(args.n_cells, len(popdf)))
            sweepdf = sweep_cells(
                sampledf,
                parameter_range=pr,
                **sweepkwargs
            )
            for i, (_, sweep) in enumerate(sweepdf.groupby(sweepdf.index)):
                sample_bitfreqs[(i, time)] = estimate_bitfreqs(sweep, dm)
                pop_repeated_bitfreqs[(i, time)] = pop_bitfreqs

            if time <= ancdf["time"].iloc[-1]:
                anc_sweepdf = sweep_cell(ancdf.loc[time].iloc[0], pr, **sweepkwargs)
                anc_bitfreqs[time] = estimate_bitfreqs(anc_sweepdf, dm)
                anc_pop_bitfreqs[time] = pop_bitfreqs

            pbar.update()

        popdatadf = make_datadf(
            sample_bitfreqs,
            sample_bitfreqs if args.self_gamma else pop_repeated_bitfreqs,
            jmed,
            jalpha,
            jweights
        ).rename_axis(["i", "time"])
        ancdatadf = make_datadf(
            anc_bitfreqs,
            anc_bitfreqs if args.self_gamma else anc_pop_bitfreqs,
            jmed,
            jalpha,
            jweights
        ).rename_axis("time")
        plot_strategy(popdatadf, ancdatadf, args.outfile, args.marker_size, args.line_width)

        logger.info("Finished")

    parser = argparse.ArgumentParser(
        description="Plot data about the frequency of strategies in a simulation."
    )
    parser.add_argument("datadir",
                        help="Directory containing the cell CSV files")
    parser.add_argument("outfile",
                        help="Output HTML or SVG file")
    parser.add_argument("-t",
                        "--n-times",
                        help="How many time points will be sampled (default: %(default)s)",
                        default=10,
                        type=int)
    parser.add_argument("-n",
                        "--n-cells",
                        help="How many cells to sample for the population strategy plot (default: %(default)s)",
                        default=100,
                        type=int)
    parser.add_argument("-g",
                        "--gamma-thresh",
                        help="Gamma threshold for classifying a strategy as propagule-former "
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
    parser.add_argument("--marker-size",
                        help="Size of the markers in the ancestry plot "
                             "(default: %(default)s)",
                        default=8,
                        type=float)
    parser.add_argument("--line-width",
                        help="Width of the lines in the ancestry plot "
                             "(default: %(default)s)",
                        default=5,
                        type=float)
    parser.add_argument("--stop-anc",
                        help="Cuts the range of the plot at the MRCA of the population (will still sample '-t' points)",
                        action="store_true")
    parser.add_argument("--Jweights",
                        help="Comma-separated list of weights used to calculate gamma between bitstrings "
                             "(default: %(default)s)",
                        default="1,2,3,4,5,6,7,8")
    parser.add_argument("--n-processes",
                        help="How many processes are used to run the analysis"
                             "(default: %(default)s)",
                        default=1,
                        type=int)
    parser.add_argument("--self-gamma",
                        help="If set, uses the ancestor gamma with itself for the plot",
                        action="store_true")
    add_inestimable_sweep_args(parser)
    parser.set_defaults(run=run)
    return parser


def plot_strategy(popdatadf, ancdatadf, outfile, anc_marker_size=15, anc_line_width=10):
    fig = go.Figure()
    strat_colors = cl.Palette.load("strats")
    popdatadf = popdatadf.sort_values(["strat", "div_gamma", "time"], ascending=[False, True, True])
    popdatadf = popdatadf.assign(
        strat_color=[get_strat_color(
            strat,
            strat_colors,
            div_gamma
        ) for strat, div_gamma in zip(popdatadf["strat"], popdatadf["div_gamma"])]
    )
    ancdatadf = ancdatadf.assign(
        strat_color=[get_strat_color(
            strat,
            strat_colors,
            div_gamma
        ) for strat, div_gamma in zip(ancdatadf["strat"], ancdatadf["div_gamma"])]
    )
    strat_count = popdatadf.groupby(level=1, as_index=False, sort=False)["strat_color"].value_counts(
        normalize=True,
        sort=False
    )
    strat_mat = pd.pivot_table(
        strat_count,
        values="proportion",
        index="strat_color",
        columns="time",
        observed=False,
        sort=False
    )
    for strat_color, row in strat_mat.iterrows():
        fig.add_trace(go.Scatter(
            x=row.index,
            y=row,
            fillcolor=strat_color,
            stackgroup="strat_color",
            name=strat_color,
            mode="lines",
            line_width=0.5,
            line_color=strat_color,
            showlegend=False
        ))

    times = ancdatadf.index.get_level_values("time")
    strat_it = zip(
        times[:-1],
        times[1:],
        ancdatadf["strat"][:-1],
        ancdatadf["strat"][1:],
        ancdatadf["strat_color"][:-1],
        ancdatadf["strat_color"][1:]
    )
    correction = 0.001 * (times[-1] - times[0])
    for prev_time, time, prev_strat, strat, prev_strat_color, strat_color in strat_it:
        if prev_strat == strat:
            fig.add_trace(go.Scatter(
                # The correction term is to fix rendering bugs with the svg
                x=[prev_time, time + correction],
                y=[-0.15, -0.15],
                mode="lines",
                line_color=strat_color,
                line_width=anc_line_width,
                showlegend=False
            ))
        else:
            fig.add_traces([
                go.Scatter(
                    x=[prev_time],
                    y=[-0.15],
                    mode="markers",
                    marker_color=prev_strat_color,
                    marker_size=anc_marker_size,
                    opacity=int(anc_marker_size != 0),
                    showlegend=False
                ),
                go.Scatter(
                    x=[time],
                    y=[-0.15],
                    mode="markers",
                    marker_color=strat_color,
                    marker_size=anc_marker_size,
                    opacity=int(anc_marker_size != 0),
                    showlegend=False
                )
            ])
    # Fake spine line
    fig.add_trace(go.Scatter(
        x=popdatadf.index.get_level_values("time").sort_values()[[0, -1]],
        y=[-0.3, -0.3],
        mode="lines",
        line_color="black",
        line_width=2,
        showlegend=False
    ))
    fig.add_traces([
        go.Scatter(
            x=[times[0]],
            y=[-0.15],
            mode="markers",
            marker_color=ancdatadf["strat_color"].iloc[0],
            marker_size=anc_marker_size,
            opacity=int(anc_marker_size != 0),
            showlegend=False
        ),
        go.Scatter(
            x=[times[-1]],
            y=[-0.15],
            mode="markers",
            marker_color=ancdatadf["strat_color"].iloc[-1],
            marker_size=anc_marker_size,
            opacity=int(anc_marker_size != 0),
            showlegend=False
        )
    ])

    fig.update_xaxes(
        title="time",
        showline=False,
        showgrid=False,
        zeroline=False,
        ticks="outside",
        anchor="free",
        position=0.01,
    )
    fig.update_yaxes(
        title="freq. of strategy",
        range=[-0.3, 1],
        showline=False,
        showgrid=False,
        zeroline=False,
        ticks="",
        anchor="free",
        position=0.06,
    )
    fig.update_layout(template="plotly_white", width=600, height=300)
    write_plot(fig, outfile)