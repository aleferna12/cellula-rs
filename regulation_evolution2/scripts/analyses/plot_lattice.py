import argparse
import logging
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from colorir import *
from scripts.fileio import *

logger = logging.getLogger(__name__)


def get_parser():
    def run(args):
        latdf = parse_lattice(args.inputfile)
        plot_lattice(latdf, args.outputfile)
        logger.info("Finished")

    parser = argparse.ArgumentParser(
        description="Plot sigma lattice."
    )
    parser.add_argument("inputfile", help="Lattice file")
    parser.add_argument("outputfile",
                        nargs='?',
                        default=None,
                        help="HTML output file (omit to show the figure in your browser instead of saving it)")
    parser.set_defaults(run=run)
    return parser


def plot_lattice(latdf: pd.DataFrame, outputfile: str = None):
    n_cells = len(np.unique(latdf.values))
    colorscale = PolarGrad(StackPalette.load("carnival"), hue_lerp="longest").to_plotly_colorscale(n_cells)
    colorscale[0] = (0.0, "#ffffff")
    fig = go.Figure(
        data=go.Heatmap(
            z=latdf.values[::-1],
            x=latdf.columns,
            y=latdf.index,
            colorscale=colorscale,
            showscale=False
        )
    )
    fig.update_layout(xaxis_constrain="domain", yaxis_scaleanchor="x", yaxis_scaleratio=1)

    if outputfile is None:
        fig.show()
    else:
        fig.write_html(outputfile)
