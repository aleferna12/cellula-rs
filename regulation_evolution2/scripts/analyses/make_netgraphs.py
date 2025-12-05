import argparse
import logging
import pickle
import networkx as nx
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from enlighten import Counter
from scripts.fileio import parse_cell_data, gene_attrs
from scripts.sweep import parse_sweep_args, add_sweep_args

logger = logging.getLogger(__name__)


def get_parser():
    def run(args):
        celldf = parse_cell_data(args.datafile)
        sweepkwargs, _ = parse_sweep_args(args)
        netgraphs = make_netgraphs(celldf, prune_alg=args.prune_level, **sweepkwargs)
        with open(args.outputfile, "wb") as file:
            pickle.dump(netgraphs, file)

    parser = ArgumentParser(description="Create network files to be plotted with 'plot_netgraph'.")
    parser.add_argument("datafile", help="Input CSV file")
    parser.add_argument("outputfile", help="Output pickle file containing the network graphs")
    add_prune_args(parser)
    parser.set_defaults(run=run)
    return parser


def read_netgraphs_file(filepath):
    with open(filepath, "rb") as file:
        return pickle.load(file)


def make_netgraphs(celldf, *args, **kwargs):
    logger.info("Creating graphs from the cells' GRNs")

    netgraphs = []
    pbar = Counter(total=len(celldf.index), desc="Networks parsed")
    for index in celldf.index:
        netgraphs.append(make_netgraph(celldf.loc[index], *args, **kwargs))
        pbar.update()

    return netgraphs


def exec_prune(netgraph):
    remove = []
    for node, pruner in netgraph.nodes(data="prune_reason"):
        if pruner == "off":
            remove.append(node)
        elif pruner == "on":
            es = netgraph.out_edges(node, data="weight")
            for e in es:
                netgraph.nodes[e[1]]["threshold"] -= e[2]
            remove.append(node)
    for node in remove:
        ntype = netgraph.nodes[node]["type"]
        netgraph.graph[ntype + "nr"] -= 1
    netgraph.remove_nodes_from(remove)


def safe_prune(netgraph):
    nx.set_node_attributes(netgraph, None, "prune_reason")
    for node, thres in netgraph.nodes(data="threshold"):
        in_es = netgraph.in_edges(node, data="weight")
        if not in_es:
            continue

        in_ws = np.array([e[2] for e in in_es])
        signs = np.sign(in_ws)
        if np.all(signs == signs[0]):
            if np.sum(in_ws) < thres:
                netgraph.nodes[node]["prune_reason"] = "off"
            else:
                netgraph.nodes[node]["prune_reason"] = "on"
    exec_prune(netgraph)


def sweep_prune(netgraph, sweepdf):
    nx.set_node_attributes(netgraph, None, "prune_reason")
    regs = np.array([np.fromstring(bit, sep=" ") for bit in sweepdf["reg_states"]], dtype=int)

    on = np.logical_and.reduce(regs, axis=0)
    off = ~np.logical_or.reduce(regs, axis=0)
    innr = netgraph.graph["innr"]
    for i, (noff, non) in enumerate(zip(off, on)):
        if noff:
            netgraph.nodes[i + innr]["prune_reason"] = "off"
        elif non:
            netgraph.nodes[i + innr]["prune_reason"] = "on"
    exec_prune(netgraph)


def count_mut(ancss: pd.Series, childss: pd.Series):
    muts = {}
    for key in gene_attrs:
        anc_a = np.fromstring(ancss.loc[key], sep=" ", dtype=float)
        ch_a = np.fromstring(childss.loc[key], sep=" ", dtype=float)
        muts[key] = np.sum(anc_a != ch_a)
    return muts


def make_netgraph(cellss: pd.Series, prune_alg=None, sweepdf=None):
    """Make directed network graph from cell series."""
    genes = {}
    for key in gene_attrs:
        genes[key] = np.fromstring(cellss[key], sep=' ', dtype=float)

    # Not doing this conversion trips up networkx bc of numpy integer and float types
    sigma, innr, regnr, outnr = cellss.loc[["sigma", "innr", "regnr", "outnr"]].astype(int)
    netgraph = nx.DiGraph(
        sigma=sigma,
        innr=innr,
        regnr=regnr,
        outnr=outnr,
    )

    # ntype is used to plot and must be numerical, but is the same as type
    for i in range(innr):
        netgraph.add_node(i, type="in", ntype=0, scale=genes["in_scale_list"][i])

    for i in range(regnr):
        reg_id = innr + i
        netgraph.add_node(reg_id, type="reg", ntype=1, threshold=genes["reg_threshold_list"][i])
        for j in range(innr):
            w_index = i * innr + j
            netgraph.add_edge(j, reg_id, weight=genes["reg_w_innode_list"][w_index])
        for j in range(regnr):
            reg_id2 = innr + j
            w_index = i * regnr + j
            netgraph.add_edge(reg_id2, reg_id, weight=genes["reg_w_regnode_list"][w_index])

    for i in range(outnr):
        out_id = innr + regnr + i
        netgraph.add_node(out_id, type="out", ntype=2, threshold=genes["out_threshold_list"][i])
        for j in range(regnr):
            reg_id = innr + j
            w_index = i * regnr + j
            netgraph.add_edge(reg_id, out_id, weight=genes["out_w_regnode_list"][w_index])

    # Pruning
    if prune_alg in ["safe", "full"]:
        safe_prune(netgraph)
    if prune_alg in ["sweep", "full"]:
        sweep_prune(
            netgraph,
            sweepdf
        )
    return netgraph


def network_to_cellss(netgraph: nx.DiGraph):
    # Reorder edges to match original format
    edges = sorted(netgraph.edges(data="weight"), key=lambda e: (e[1], e[0]))
    genome_data = {
        "innr": netgraph.graph["innr"],
        "regnr": netgraph.graph["regnr"],
        "outnr": netgraph.graph["outnr"],
        "in_scale_list": ' '.join(str(x) for x in nx.get_node_attributes(netgraph, "scale").values()),
        "reg_threshold_list": ' '.join(
            str(data["threshold"]) for _, data in netgraph.nodes(data=True)
            if data["type"] == "reg"
        ),
        "reg_w_innode_list": ' '.join(
            str(edge[2]) for edge in edges
            if netgraph.nodes[edge[1]]["type"] == "reg"
            and netgraph.nodes[edge[0]]["type"] == "in"
        ),
        "reg_w_regnode_list": ' '.join(
            str(edge[2]) for edge in edges
            if netgraph.nodes[edge[0]]["type"] == netgraph.nodes[edge[1]]["type"]
        ),
        "out_threshold_list": ' '.join(
            str(data["threshold"]) for _, data in netgraph.nodes(data=True)
            if data["type"] == "out"
        ),
        "out_w_regnode_list": ' '.join(
            str(edge[2]) for edge in edges
            if netgraph.nodes[edge[1]]["type"] == "out"
        ),
    }
    return pd.Series(genome_data)


def add_prune_args(parser: argparse.ArgumentParser):
    prune = parser.add_argument_group(
        "pruning arguments",
        description="Arguments specifying how to prune the networks. If 'prune-level' is 'sweep' or 'full', "
                    "sweep arguments must also be specified."
    )
    prune.add_argument("--prune-level",
                       help="Algorithm used for pruning (default: %(default)s)",
                       choices=["none", "safe", "sweep", "full"],
                       default="none")
    add_sweep_args(parser)
