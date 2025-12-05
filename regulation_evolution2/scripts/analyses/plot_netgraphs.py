import argparse
import logging
import networkx as nx
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import graphviz as gv
import pandas as pd
from enlighten import Counter
from colorir import *
from scripts.analyses.make_netgraphs import *
from scripts.fileio import parse_cell_data
from scripts.config import CS

logger = logging.getLogger(__name__)
_node_cs = {
    "in": (CS.redlipstick, blend(CS.redblush, CS.white, 0.7)),
    "reg": (CS.deeplagoon, blend(CS.shallowsea, CS.white, 0.7)),
    "out": (CS.matcha, blend(CS.greentea, CS.white, 0.7))
}


def get_parser():
    def run(args):
        if args.datafile[-4:].lower() == ".csv":
            celldf = parse_cell_data(args.datafile)
            if args.sigmas is not None:
                celldf = celldf[celldf["sigma"].isin(np.array(args.sigmas, dtype=int))]
            sweepkwargs, _ = parse_sweep_args(args)
            netgraphs = make_netgraphs(celldf, prune_alg=args.prune_level, **sweepkwargs)
        else:
            netgraphs = read_netgraphs_file(args.datafile)
        pbar = Counter(total=len(netgraphs), desc="Networks plotted")
        for netgraph in netgraphs:
            plot_netgraph(netgraph,
                          f"{args.outputdir}/{netgraph.graph['sigma']}.svg",
                          args.color_prune,
                          args.reg_edges)
            pbar.update()
        logger.info("Finished")

    parser = argparse.ArgumentParser(
        description="Plot GRNs as graphs."
    )
    parser.add_argument("datafile",
                        help="Either a pickle or a CSV file to read the networks from")
    parser.add_argument("outputdir", help="directory where to save the SVGs")
    parser.add_argument("-p",
                        "--color-prune",
                        help="Visual \"prune\" level, the higher the less edges will be shown "
                             "(default: %(default)s)",
                        default=1,
                        type=int)
    parser.add_argument("-r",
                        "--reg-edges",
                        help="Whether to plot edges between regulation nodes",
                        action="store_true")
    parser.add_argument("-b",
                        "--backend",
                        help="Whether to plot the network with matplotlib or "
                             "graphviz (default: %(default)s)",
                        choices=["matplotlib", "graphviz"],
                        default="matplotlib")
    parser.add_argument("-s",
                        "--sigmas",
                        help="List of space-delimited sigmas to plot",
                        default=None,
                        nargs='*')
    add_prune_args(parser)
    parser.set_defaults(run=run)
    return parser


def plot_netgraph(netgraph, outfile, color_prune=1, reg_edges=True, backend="matplotlib"):
    netgraph = netgraph.copy()

    for node, ntype in netgraph.nodes(data="type"):
        netgraph.nodes[node]["node_color"] = _node_cs[ntype][1]
        netgraph.nodes[node]["edge_color"] = _node_cs[ntype][0]
        if ntype == "in":
            continue

        in_ws = np.array([e[2] for e in netgraph.in_edges(node, data="weight")])
        netgraph.nodes[node]["weight_range"] = np.sum(in_ws, where=in_ws > 0) - np.sum(in_ws, where=in_ws < 0)

    weights = netgraph.edges(data="weight")
    relativ_weights = []
    for edge in weights:
        relativ_weights.append(edge[2] / netgraph.nodes[edge[1]]["weight_range"])
    grad = Grad(
        [CS.redlipstick]
        + [sRGB(0, 0, 0, 0)] * color_prune
        + [CS.matcha],
        domain=[-1, 1],
        color_format=MATPLOTLIB_COLOR_FORMAT
    )

    edge_colors = [grad(w) for w in relativ_weights]
    nx.set_edge_attributes(netgraph, dict(zip(netgraph.edges, edge_colors)), "color")

    if backend == "graphviz":
        return draw_netgraph_graphviz(netgraph, outfile, reg_edges)
    else:
        return draw_netgraph_plt(netgraph, outfile, reg_edges)


def multipartite_layout(netgraph):
    xs = list(nx.get_node_attributes(netgraph, "ntype").values())
    ys = np.concatenate([np.linspace(-1, 1, n + 2)[1:-1]
                         for n in [netgraph.graph["innr"],
                                   netgraph.graph["regnr"],
                                   netgraph.graph["outnr"]]])
    return dict(zip(netgraph.nodes, np.stack([xs, ys], axis=1)))


def draw_netgraph_plt(netgraph, outfile, reg_edges):
    fig = plt.figure(figsize=(6, 8))
    pos = multipartite_layout(netgraph)
    nx.draw_networkx_labels(netgraph, pos, font_size=8)
    nx.draw_networkx_nodes(netgraph,
                           pos,
                           node_size=200,
                           linewidths=1.5,
                           node_color=list(nx.get_node_attributes(netgraph, "node_color").values()),
                           edgecolors=list(nx.get_node_attributes(netgraph, "edge_color").values()))

    es = {}
    for e in netgraph.edges(data="color"):
        # Both are regnodes
        if netgraph.nodes[e[0]]["type"] == netgraph.nodes[e[1]]["type"]:
            if not reg_edges:
                continue
            if e[0] == e[1]:
                conn = "arc3"
            else:
                conn = "arc3,rad=0.4"
        else:
            conn = "arc3"
        et = es.setdefault(conn, ([], []))
        et[0].append(e)
        et[1].append(e[2])

    for conn, edges in es.items():
        nx.draw_networkx_edges(netgraph, pos, edgelist=edges[0], edge_color=edges[1], connectionstyle=conn)
    plt.savefig(outfile)
    plt.close()
    return fig


def draw_netgraph_graphviz(netgraph, outfile, reg_edges):
    g = gv.Digraph()
    g.attr(ranksep="5", rankdir="LR", nodesep="0.2")
    g.attr('node', shape='circle', style='filled', penwidth="2")
    for node, data in netgraph.nodes(data=True):
        g.node(f"{node:02}", fillcolor=data["node_color"], color=data["edge_color"])
    g.attr("edge", penwidth="2")
    for edge in netgraph.edges(data="color"):
        con = "true"
        if netgraph.nodes[edge[0]]["type"] == netgraph.nodes[edge[1]]["type"]:
            if not reg_edges:
                continue
            con = "false"
        g.edge(
            f"{edge[0]:02}",
            f"{edge[1]:02}",
            constraint=con,
            color=edge[2].hex(include_a=True, tail_a=True)
        )
    g.render(outfile=outfile, cleanup=True)
    return g
