import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scripts.calculate_adh import *

POP = 400000
JMED = 24
JALPHA = 12
MULT = np.arange(1, 9)
MULT = np.array(MULT)


def main():
    jkeys = make_bitstrings(POP, len(MULT))
    jlocks = make_bitstrings(POP, len(MULT))
    ccJs = []
    gammas = []
    it = enumerate(zip(jkeys, jlocks))
    adhtable = contact_energy_table(MULT)
    for i, (jkey1, jlock1) in it:
        for j, (jkey2, jlock2) in it:
            if i != j:
                ccJ = cell_contact_energy(jkey1, jlock1, jkey2, jlock2, adhtable)
                ccJs.append(ccJ)
                gammas.append(calculate_gamma(JMED, JALPHA, ccJ))

    fig = make_subplots(2, 1)
    values, counts = np.unique(np.round(ccJs, 8), return_counts=True)
    fig.add_trace(go.Scatter(x=values, y=counts / np.sum(counts)), row=1, col=1)
    values, counts = np.unique(np.round(gammas, 8), return_counts=True)
    fig.add_trace(go.Scatter(x=values, y=counts / np.sum(counts)), row=2, col=1)
    fig.update_layout(xaxis_title="Jcc",
                      yaxis_title="relative freq.",
                      title="Weights: " + " ".join(str(x) for x in MULT))
    fig.update_yaxes(range=[0, 0.25])
    fig.update_xaxes(range=[0, 50], row=1, col=1)
    fig.update_xaxes(range=[-25, 25], row=2, col=1, title="gamma")
    # fig.write_html(Path("~/Desktop/Jcc.html").expanduser())
    fig.show()


def make_bitstrings(n, length):
    rng = np.random.default_rng()
    return rng.integers(0, 2 ** length, n)


if __name__ == "__main__":
    main()
