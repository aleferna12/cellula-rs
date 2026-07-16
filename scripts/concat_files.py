#!/usr/bin/env python3
import sys
from pathlib import Path

import polars as pl


def concat_dir(dir_path: Path) -> pl.DataFrame:
    if not any(dir_path.glob("*.parquet")):
        raise ValueError(f"No .parquet files found in {dir_path}")

    df = pl.read_parquet(dir_path / "*.parquet", include_file_paths="file_path")
    df = df.with_columns(
        wtime=pl.col("file_path")
        .str.split("/")
        .list.get(-1)
        .str.replace(".parquet", "", literal=True)
        .cast(pl.UInt32)
    )
    return df.drop("file_path")


def main():
    if len(sys.argv) != 2:
        print("Usage: python concat_parquet.py <root_dir>", file=sys.stderr)
        sys.exit(1)

    root = Path(sys.argv[1])
    cells_dir = root / "cells"
    lattices_dir = root / "lattices"

    if not cells_dir.is_dir() or not lattices_dir.is_dir():
        print(
            f"Error: both 'cells' and 'lattices' subdirectories must exist under {root}",
            file=sys.stderr,
        )
        sys.exit(1)

    for name, dir_path in (("cells", cells_dir), ("lattices", lattices_dir)):
        df = concat_dir(dir_path)
        out_path = root / f"{name}.parquet"
        df.write_parquet(out_path)
        print(f"Wrote {out_path} ({df.height} rows)")


if __name__ == "__main__":
    main()
