import polars as pl
from pathlib import Path


def read_celldfs(top_path, levels=["replica", "energy"]):
    df = pl.read_parquet(Path(top_path) / "**" / "cells" / "*.parquet", include_file_paths="file_path")
    df = df.with_columns(path_list=pl.col("file_path").str.split("/"))

    if levels:
        ldict = {lv: pl.col("path_list").list.get(-i - 3) for i, lv in enumerate(levels)}
        df = df.with_columns(
            wtime=pl.col("path_list").list.get(-1).str.replace(".parquet", "").cast(pl.UInt32),
            **ldict
        )
    return df.select(pl.exclude(["path_list", "file_path"]))


def save_plot(fig, name):
    fig.write_html(f"{name}.html")
    fig.write_image(f"{name}.svg")
    return fig