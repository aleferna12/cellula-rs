import random

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scripts.calculate_adh import *
from scripts.fileio import write_plot


def main():
    # plot_perfect_bit(14, 7, np.arange(1, 7))
    plot_selection(5, 14, 7, np.arange(1, 7))


def plot_selection(gamma_thresh, jmed, jalpha, jweights):
    ref_bits = []
    adht = contact_energy_table(jweights)
    for i in range(2 ** len(jweights)):
        for j in range(2 ** len(jweights)):
            jcc = cell_contact_energy(i, j, i, j, adht)
            gamma = calculate_gamma(jmed, jalpha, jcc)
            if gamma >= gamma_thresh:
                ref_bits.append([gamma, i, j])

    dfs = []
    for gamma, key, lock in random.sample(ref_bits, min(100, len(ref_bits))):
        adhdf = make_adhdf(key, lock, jmed, jalpha, jweights, 10)
        adhdf["bit"] = f"{key} {lock} - " + adhdf["bit"]
        dfs.append(adhdf)
    df = pd.concat(dfs)
    fig = px.scatter(df.sort_values(["neg", "type", "gamma"]), x="bit", y="gamma", color="type", symbol="neg")
    fig.show()


def plot_perfect_bit(jmed, jalpha, jweights):
    rkey = np.random.rand(len(jweights)) < 0.5
    rlock = ~rkey
    df = make_adhdf(bitarray_to_dec(rkey), bitarray_to_dec(rlock), jmed, jalpha, jweights)
    fig = px.scatter(df.sort_values(["type", "gamma"]), x="bit", y="gamma", color="type", symbol="neg")
    fig.show()


def make_adhdf(rkey, rlock, jmed, jalpha, jweights, n=None):
    adht = contact_energy_table(jweights)
    gammas = {"bit": [], "gamma": [], "type": [], "neg": []}
    if n is None:
        iti = range(2 ** len(jweights))
        itj = range(2 ** len(jweights))
    else:
        iti = random.sample(range(2 ** len(jweights)), n)
        itj = random.sample(range(2 ** len(jweights)), n)
    for i in iti:
        for j in itj:
            rjcc = cell_contact_energy(rkey, rlock, i, j, adht)
            sjcc = cell_contact_energy(i, j, i, j, adht)
            rgamma = calculate_gamma(jmed, jalpha, rjcc)
            sgamma = calculate_gamma(jmed, jalpha, sjcc)
            for gamma in [rgamma, sgamma, (rgamma + sgamma) / 2]:
                gammas["bit"].append(str(i) + ' ' + str(j))
                gammas["gamma"].append(gamma)
            gammas["type"].extend(["ref", "self", "avg"])
            gammas["neg"].extend([rgamma < 0 and sgamma < 0] * 3)

    return pd.DataFrame.from_dict(gammas)


if __name__ == "__main__":
    main()

