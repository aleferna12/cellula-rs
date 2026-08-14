import polars as pl
from pathlib import Path


def read_celldfs(top_path, levels=["replica", "energy"], low_memory=False, concated=False, scan=False):
    top_path = Path(top_path)
    if concated:
        search_path = top_path / "**" / "cells.parquet"
    else:
        search_path = top_path / "**" / "cells" / "*.parquet"
    method = pl.scan_parquet if scan else pl.read_parquet
    df = pl.read_parquet(search_path, include_file_paths="file_path", low_memory=low_memory)
    df = df.with_columns(path_list=pl.col("file_path").str.split("/"))

    if levels:
        ldict = {lv: pl.col("path_list").list.get(-i - 3 + concated) for i, lv in enumerate(levels)}
        if not concated:
            df = df.with_columns(
                wtime=pl.col("path_list").list.get(-1).str.replace(".parquet", "").cast(pl.UInt32),
            )
        df = df.with_columns(
            **ldict
        )
    return df.select(pl.exclude(["path_list", "file_path"]))


def save_plot(fig, name):
    fig.write_html(f"{name}.html")
    fig.write_image(f"{name}.svg")
    return fig