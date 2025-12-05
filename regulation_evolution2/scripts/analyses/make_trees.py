"""Makes newick trees from cell ancestry data.

Be mindful of the fact that not all generations are captured, since data is only saved every X
MCSs. This means that some branches will be collapsed into multifurcations, resulting in some loss
of phylogenetic information. This can be improved (with some amount of work) by saving a vector of
ancestors rather than a single one.
"""
import argparse
import logging
import pandas as pd
from pathlib import Path
from enlighten import Counter
from ete3 import Tree
from scripts.fileio import *
logger = logging.getLogger(__name__)


def get_parser():
    def run(args):
        t_filter = build_time_filter(get_time_points(args.datadir), end=args.last_time)
        celldf = parse_cell_data(args.datadir, time_filter=t_filter)
        trees = make_trees(celldf,
                           extinct=args.extinct,
                           names=not args.no_names,
                           single_lineage=args.single_lineage,
                           stop_mrca=args.stop_mrca,
                           collapse_branches=not args.no_collapse)

        longest = get_longest_trees(trees)
        logger.info(f"Longest lineages survived for {longest[0][1]} MCSs)")

        logger.info(f"Writing trees to '{args.outdir}'")
        outdir = Path(args.outdir)
        for i, (tree, _) in enumerate(longest):
            tree_name = f"tree{i}.newick" if len(trees) > 1 else args.outfile
            tree.write(
                outfile=str(outdir / tree_name),
                format=args.format,
                features=None if args.no_nhx else ["time"],
                format_root_node=not args.no_root
            )

        logger.info("Finished")

    parser = argparse.ArgumentParser(
        description="Create tree files to be plotted with 'plot_tree'."
    )
    parser.add_argument("datadir", help="Directory containing the cell data as CSV files")
    parser.add_argument("outdir", help="Directory where the trees will be created")
    parser.add_argument("-e", 
                        "--extinct", 
                        action="store_true",
                        help="Include extinct lineages")
    parser.add_argument("-n", 
                        "--no-names",
                        action="store_true",
                        help="Do not include the sigmas of cells as branch names")
    parser.add_argument(
        "-s",
        "--single-lineage",
        action="store_true",
        help="Allow trees that consist of a single lineage to be created (may cause problems with "
             "visualization in some programs)"
    )
    parser.add_argument("-x",
                        "--no-nhx",
                        action="store_true",
                        help="Do not include NHX attributes in the trees (such as time)")
    parser.add_argument("-m",
                        "--stop-mrca",
                        action="store_true",
                        help="Stop computing when the MRCA of the living cells is found")
    parser.add_argument(
        "-c",
        "--no-collapse",
        action="store_true",
        help="Do not collapse lineages with a single descendent into a single branch"
    )
    parser.add_argument("-r",
                        "--no-root",
                        action="store_true",
                        help="Do not include root attributes and information for compatibility")
    parser.add_argument("-f",
                        "--format",
                        default=5,
                        type=int,
                        help="Tree output format (following etetoolkit conventions) (default: %(default)s)")
    parser.add_argument("-t",
                        "--last-time",
                        default=-1,
                        type=int,
                        help="Last time point to look at when making the trees (i.e. the present)")
    parser.add_argument("-o",
                        "--outfile",
                        help="Specifies the name of the file used when the analysis results in a "
                             "single tree (otherwise they will be numbered) (default: %(default)s)",
                        default="tree.newick")
    parser.set_defaults(run=run)
    return parser


def get_longest_trees(trees: list[Tree]):
    tree_times = [(tree, int(tree.get_farthest_leaf()[0].time) - int(tree.time)) for tree in trees]
    return sorted(tree_times, key=lambda x: x[1], reverse=True)


# Remember that although the sigmas are added as names, the same id can refer to DIFFERENT cells
def make_trees(celldf: pd.DataFrame,
               extinct=False,
               names=True,
               single_lineage=True,
               stop_mrca=False,
               collapse_branches=True) -> list[Tree]:
    logger.info("Building trees from cell ancestry data")
    # This dict holds the descendents of the current cells
    prev_anc_child = {}
    # Order by time backwards
    # This is needed to deal with dead-ends in a nice way and to terminate early if MRCA is found
    timesteps = celldf["time"].sort_values(ascending=False).unique()
    if timesteps.size == 0:
        raise ValueError("'celldf' argument must have more than one timestep for ancestry tracing")

    pbar = Counter(total=len(timesteps), desc="Time-steps")
    for time in timesteps:
        # This dict holds the ancestors of the current cells
        next_anc_child = {}
        tdf = celldf.loc[time]
        for sigma, anc_sigma in zip(tdf["sigma"], tdf["ancestor"]):
            # Keep dead ends from being added to the trees
            if not extinct and prev_anc_child and sigma not in prev_anc_child:
                continue

            node = _build_node(
                sigma,
                prev_anc_child.get(sigma, []),
                time,
                names,
                collapse_branches
            )
            if anc_sigma not in next_anc_child:
                next_anc_child[anc_sigma] = []
            next_anc_child[anc_sigma].append(node)

        found_mrca = len(prev_anc_child) == 1
        prev_anc_child = dict(next_anc_child)
        if stop_mrca and found_mrca:
            logger.info("Found MRCA, stopping analysis")
            pbar.close()
            break
        pbar.update()

    # Construct root trees
    trees = []
    for root_sigma, roots in prev_anc_child.items():
        if len(roots) > 1:
            raise Exception(f"found multiple roots for tree {root_sigma}: {roots}. Something is (terribly) wrong!")
        root = roots[0]
        # Exclude trees that only have a single leaf due to being part of an exclusively
        # migrating lineage or because the lineage only has one survivor and extinct=False
        # These single-leaf lineages can cause problems in some newick parsers
        if single_lineage or len(root.get_leaves()) > 1:
            trees.append(root)
    return trees


def _build_node(sigma, children, time, names, collapse_branches):
    if collapse_branches and len(children) == 1:
        children[0].dist += 1
        return children[0]
    node = Tree(name=str(sigma) if names else "")
    node.add_feature("time", time)
    for child in children:
        node.add_child(child)
    return node
