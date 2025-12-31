from pathlib import Path
from multiprocessing import get_context
import polars as pl


def parse_cells(sim_path):
    dfs = []
    for file in Path(sim_path / "cells").iterdir():
        data = pl.read_parquet(file)
        time = int(file.name.rstrip(".parquet"))
        data = data.with_columns(time=time)
        dfs.append(data)
    return pl.concat(dfs)


def parse_cells_multiple(sim_paths: list, n_workers: int):
    pool = get_context("spawn").Pool(n_workers)
    results = []
    for sim_path in sim_paths:
        res = pool.apply_async(parse_cells, (sim_path,))
        results.append(res)
    pool.close()
    pool.join()
    return [r.get() for r in results]