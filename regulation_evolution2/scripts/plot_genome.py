from pathlib import Path

import pandas as pd
import plotly.express as px
from colorir import StackPalette

from scripts.analyses.plot_regulation import chem_to_dist
from scripts.fileio import *


def main():
    for file in Path("test/output").glob("*.csv"):
        df = pd.read_csv(file)
        df = make_df(df)
        f = px.scatter(df,
                       x="fdist",
                       y="food",
                       color="tau_str",
                       symbol="bitstring",
                       color_discrete_map=dict(zip(["1", "2"], StackPalette.load("plotly")[0:2][::-1])))
        dr = df["fdist"].max() - df["fdist"].min()
        fr = df["food"].max() - df["food"].min()
        f.update_layout(xaxis_constrain="domain", yaxis_scaleanchor="x", yaxis_scaleratio=dr / fr)  # Make square
        write_plot(f, str(file).replace(".csv", ".html"))


def make_df(df: pd.DataFrame):
    df = df.copy()
    df["fdist"] = [chem_to_dist(x, 2000, 5) for x in df["chem"]]
    df["bitstring"] = df["jkey_dec"].astype(str) + "-" + df["jlock_dec"].astype(str)
    df["food"] = df["food"] * 20 * 1020 / 100
    df["tau_str"] = df["tau"].astype(str)
    return df[(df["fdist"] > 0) & (df["food"] < 3000)]


if __name__ == "__main__":
    main()

