import polars as pl
from ete3 import Tree


def matrix(celldf, season_duration):
    wtimes = celldf["wtime"].sort().diff()
    wtime = wtimes.filter(wtimes != 0).mean()
    ancs = (
        celldf
            .filter(
                (pl.col("wtime") - wtime) % season_duration == 0  # Only changes on time after season change
            )
            .select(["index", "ancestor", "wtime"])
            .sort("wtime")
    )
    ancdf = ancs.filter(wtime=ancs["wtime"].first()).select("index")
    for (time,), group in ancs.group_by("wtime", maintain_order=True):
        if time == ancs["wtime"].first():
            continue
        ancdf = (
            group
                .join(
                    ancdf,
                    left_on="ancestor", 
                    right_on="index",
                    how="full",  # Change to left to exclude extinct lineages (same as dropping nulls from the matrix afterwards)
                    coalesce=False
                )
                .select(["index", "index_right", pl.col(r"^\d*$")])
                .rename({
                    "index_right": f"{time - int(season_duration)}"
                })
        )
    ancdf = ancdf.rename({"index": str(time)})
    return ancdf


def trees(phy_matrix):
    prev_time = "0"
    nodes = [None] * phy_matrix[prev_time].n_unique()
    for x in phy_matrix[prev_time].unique():
        t = Tree(name=str(x))
        t.add_feature("time", "0")
        nodes[x] = t
    for row in phy_matrix.transpose():
        row = row[::-1]
        prev_anc = row[0]
        node = nodes[prev_anc]
        for anc in row.filter(row.is_not_null())[1:]:
            anc_str = str(anc)
            found = False
            for child in node.children:
                if child.name == anc_str:
                    found = True
                    break
            if found:
                node = child
            else:
                time = str(2_000_000 + int(node.time))
                node = node.add_child(name=anc_str)
                node.add_feature("time", time)
    return nodes