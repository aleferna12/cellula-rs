import argparse
import warnings
import logging
import re
from pathlib import Path
from ete3 import Tree, TreeStyle, AttrFace
from random import shuffle, seed
from scripts.fileio import parse_cell_data
from scripts.analyses.plot_timeline import get_adhering_clusters, get_cluster_colors, CellCluster
from scripts.config import CS

logger = logging.getLogger(__name__)


def get_parser():
    def run(args):
        logger.info("Reading tree file")
        tree = Tree(args.treepath, format=5)
        # When tree is constructed by neighbour joining they are unrooted so we must estimate
        if args.reroot:
            mid = tree.get_midpoint_outgroup()
            tree.set_outgroup(mid)

        clusters = None
        if args.clusterfile is not None:
            celldf = parse_cell_data(args.clusterfile)
            clusters = get_adhering_clusters(celldf)
        plot_tree(
            tree,
            args.outfile,
            clusters,
            args.min_cluster,
            args.show_attrs,
            scale=args.tree_scale,
            cluster_dist=args.cluster_distance
        )
        logger.info("Finished")

    parser = argparse.ArgumentParser(
        description="Plot the phylogenetic relations of the simulation."
    )
    parser.add_argument("treepath", help="NEWICK file to be plotted")
    parser.add_argument("outfile", help="Output image file")
    parser.add_argument(
        "-c",
        "--clusterfile",
        help="CSV file containing the cell data from the most recent time-step present in the tree "
             "(used to draw the clusters)"
    )
    parser.add_argument("-m",
                        "--min-cluster",
                        default=2,
                        type=int,
                        help="Minimum number of cells in a cluster for it to be assigned a color "
                             "(default: %(default)s)")
    parser.add_argument("-s",
                        "--tree-scale",
                        default=100,
                        type=float,
                        help="Tree scale for plotting (default: %(default)s)")
    parser.add_argument("-d",
                        "--cluster-distance",
                        default=10,
                        type=float,
                        help="Distance for the colored face that represents clusters (default: %(default)s)")
    parser.add_argument("-a",
                        "--show-attrs",
                        action="store_true",
                        help="Shows the name and NHX attributes associated with each node (such as 'time')")
    parser.add_argument("-r",
                        "--reroot",
                        action="store_true",
                        help="Reroot the tree on a mid-point before plotting (not recommended)")
    parser.set_defaults(run=run)
    return parser


def plot_tree(tree: Tree,
              outpath,
              clusters: list[CellCluster] = None,
              min_cluster=2,
              show_attrs=False,
              scale=100,
              cluster_dist=10):
    tree = tree.copy()
    logger.info(f"Plotting tree to '{outpath}'")

    # If root dist is too short it causes problems displaying branch attributes
    if tree.dist == 1:
        tree.dist = tree.children[0].dist / 2

    leaf_color = {}
    if clusters is not None:
        colors = [str(color.hex()) for color in get_cluster_colors(clusters, min_cluster)]
        for cluster, color in zip(clusters, colors):
            for leaf in cluster:
                leaf_color[str(leaf)] = color

    last_timepoint = getattr(tree.get_farthest_leaf()[0], "time", None)
    for node in tree.traverse():
        node.img_style["size"] = 0
        node.img_style["vt_line_width"] = 50
        node.img_style["hz_line_width"] = 50
        if node.is_leaf():
            if show_attrs:
                # rapidnj adds these for no reason
                node.name = node.name.strip("\'")
                face = AttrFace("name")
                node.add_face(face, column=0)
            # Color unicellular nodes dark gray and doesn't color dead-end nodes
            if clusters is not None and getattr(node, "time", None) == last_timepoint:
                if node.name not in leaf_color:
                    raise ValueError("'clusters' doesn't contain one or more terminal leaves "
                                     "from 'tree', check if clusters and tree come from the same "
                                     f"time point: {last_timepoint}")
                node.img_style["bgcolor"] = leaf_color[node.name]
                node.img_style["hz_line_color"] = leaf_color[node.name]
                node.dist = cluster_dist
        elif show_attrs and "time" in node.features:
            face = AttrFace("time")
            node.add_face(face, column=0, position="branch-top")

    ts = TreeStyle()
    ts.root_opening_factor = 0.1
    ts.show_leaf_name = False
    ts.mode = 'c'
    ts.scale = scale

    tree.render(str(outpath), tree_style=ts, w=100, units="mm")


def figtree_nexus_str(newick, clusters, min_cluster=2):
    """Creates a color-coded nexus string that can be parsed by FigTree.

    This approach was abandoned because the clusters couldn't be distinguished in the end result.
    That said, the distances between the branches are much easier to see in FigTree thanks to the
    "radial" view (but you don't need this function for that, just open the normal tree in FigTree
    or Dendroscope).
    """
    def sub_color(match):
        number = match.group(1)
        if number not in leaf_color:
            return match.group()
        lcolor = leaf_color.get(number, CS.black)
        return f"{number}[&!color={lcolor}]:"

    colors = get_cluster_colors(clusters, min_cluster)
    leaf_color = {}
    for cluster, color in zip(clusters, colors):
        for leaf in cluster:
            leaf_color[str(leaf)] = color
    figtree_fmt = re.sub(r"(\d+):", sub_color, newick)

    return (
        "#NEXUS\n"
        "Begin Trees;\n"
        f" Tree tree1 = [&R] {figtree_fmt}\n"
        "End;\n"
    )
