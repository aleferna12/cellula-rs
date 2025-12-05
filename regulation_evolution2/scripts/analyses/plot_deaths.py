import argparse
import logging
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from colorir import *
from scripts.fileio import *

logger = logging.getLogger(__name__)
colors = px.colors.qualitative.Vivid


def get_parser():
    def run(args):
        gravedf = parse_grave_data(args.datadir)
        plot_deaths(gravedf, args.outputfile, args.bin_size)
        logger.info("Finished")

    parser = argparse.ArgumentParser(
        description="Plot information about cell death."
    )
    parser.add_argument("datadir", help="Directory containing the cell grave data frames")
    parser.add_argument("outputfile", help="SVG output file")
    parser.add_argument("-b",
                        "--bin-size",
                        help="Determine the plot bin size in MCS (default: %(default)s)",
                        type=int,
                        default=50000)
    parser.set_defaults(run=run)
    return parser


def plot_deaths(gravedf: pd.DataFrame, outputfile, bin_size):
    bins = range(gravedf["time_death"].min(), gravedf["time_death"].max() + bin_size, bin_size)
    intervals = pd.cut(gravedf["time_death"], bins=bins)
    deaths_interval = gravedf.groupby(intervals).size()
    fmt_str = "%i - %i "
    if bin_size >= 1E6:
        fmt_str += "mMCS"
        div = 1E6
    else:
        fmt_str += "kMCS"
        div = 1E3
    x = [fmt_str % (ti.left / div, ti.right / div)
         for ti in intervals.cat.categories]

    fig = make_subplots(
        3,
        1,
        shared_xaxes=True,
        subplot_titles=["Number of deaths", "Cause of death in %", "Mean age of death"]
    )

    # Plot reason of death
    reasons = gravedf.groupby(["reason", intervals]).size()
    for i, reason in enumerate(np.sort(gravedf["reason"].unique())):
        color = colors[i]
        fig.add_trace(go.Scatter(
            x=x,
            y=reasons[reason],
            stackgroup=0,
            line_color=color,
            name=reason
        ), row=1, col=1)
        fig.add_trace(go.Scatter(
            x=x,
            y=reasons[reason] / deaths_interval,
            stackgroup=1,
            line_color=color,
            showlegend=False
        ), row=2, col=1)

    # Plot average age of death
    fig.add_trace(go.Scatter(
        x=x,
        y=gravedf["time_since_birth"].groupby(intervals).mean(),
        name="age",
        showlegend=False,
    ), row=3, col=1)

    fig.update_yaxes(range=[0, 1], row=2)
    fig.update_xaxes(range=[0, len(x) - 1])

    write_plot(fig, outputfile)
