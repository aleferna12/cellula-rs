import polars as pl

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
                    "index_right": f"{time - int(season_duration + wtime)}"
                })
        )
    ancdf = ancdf.rename({"index": str(time - int(wtime))})
    return ancdf