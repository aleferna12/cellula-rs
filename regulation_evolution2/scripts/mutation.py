import tempfile
import colorir as cl
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from multiprocessing import Pool
from scripts.fileio import *
from scripts.calculate_adh import *
from scripts.sweep import mutation_wrapper, sweep_genome_wrapper
from scripts.analyses.make_netgraphs import observed_parameter_range


def main():
    datapath = "/home/aleferna/CPM/food/a_runs/21/celldata_"
    tf = build_time_filter(get_time_points(datapath), start=10000, n=500)
    celldf = parse_cell_data(datapath, n_processes=10, time_filter=tf)
    replicas = 500
    jweights = np.arange(1, 9)
    mut_rates = [0.05]
    mut_stds = [0.05]
    par_range, _ = observed_parameter_range(celldf, include_thres=0.99)
    nruns = len(mut_rates) * len(mut_stds)
    with Pool(nruns) as pool:
        arglist = []
        for mut_rate in mut_rates:
            for mut_std in mut_stds:
                arglist.append([
                    jweights,
                    celldf.iloc[0],
                    replicas,
                    mut_rate,
                    mut_std,
                    par_range["grad_conc"].min(),
                    par_range["grad_conc"].max(),
                    2.5,
                    0,
                    par_range["food"].max() / 204,
                    0.5,
                    50,
                    True
                ])
        results = pool.starmap(evolution_replica, arglist)
    datadfs = []
    tau_diff = []
    tau_diff_counts = []
    jcc_diff = []
    jcc_diff_count = []
    fig = make_subplots(2, 1)
    colors = cl.StackPalette.load("plotly", color_format=cl.WEB_COLOR_FORMAT).resize(nruns)
    for i, (args, x) in enumerate(zip(arglist, results)):
        datadfs.append(x[0])
        rate_std = f"mut_rate={args[3]} mut_std={args[4]}"
        trace_name = f"{i}: {rate_std}"
        fig.add_trace(go.Scatter(
            mode="lines",
            y=x[1],
            legendgroup=rate_std,
            name=trace_name,
            line_color=colors[i]
        ), 1, 1)
        fig.add_trace(go.Violin(
            y=x[1],
            legendgroup=rate_std,
            name=trace_name,
            line_color="black",
            fillcolor=colors[i],
            opacity=0.6
        ), 2, 1)
        tau_diff.append(x[2])
        tau_diff_counts.append(x[3])
        jcc_diff.append(x[4])
        jcc_diff_count.append(x[5])
    axes = datadfs[0]["chem"].unique(), datadfs[0]["foodp"].unique()[::-1]

    fig.show("browser")
    for data, max_val in zip([tau_diff, jcc_diff], [1, 2 * np.sum(jweights)]):
        px.imshow(
            np.array(data),
            x=axes[0],
            y=axes[1],
            facet_col=0,
            facet_col_wrap=len(mut_stds),
            animation_frame=1,
            range_color=[0, max_val],
            aspect="auto",
            origin="lower"
        ).show("browser")
    for data in [tau_diff_counts, jcc_diff_count]:
        px.imshow(
            np.array(data),
            x=axes[0],
            y=axes[1],
            facet_col=0,
            facet_col_wrap=len(mut_stds),
            aspect="auto",
            origin="lower"
        ).show("browser")
    pass


def evolution_replica(jweights, *datadf_args):
    datadf = make_datadf(*datadf_args)

    tau_matrix = tau_evolution(datadf)
    tau_diff = np.mean(tau_matrix != tau_matrix[0], axis=(1, 2), where=~np.isnan(tau_matrix))
    anc_diffs = np.where(np.isnan(tau_matrix), tau_matrix, tau_matrix != tau_matrix[0])
    anc_diff_counts = np.sum(anc_diffs, axis=0)

    jk_matrix, jl_matrix = bit_evolution(datadf)
    bit_diff = (jk_matrix != jk_matrix[0]) | (jl_matrix != jl_matrix[0])
    indexes = np.argwhere(bit_diff)
    # We can use the Jcc logic to compare the bitstrings
    # However, we are not interested in the Jcc values, we only want to quantify how much of an impact on fitness
    # a random mutation would have (which correlates with the weight of the gene being flipped)
    weighted_diff = np.zeros_like(bit_diff, float)
    adh_table = contact_energy_table(jweights)
    max_diff = cell_contact_energy(0, 0, 0, 0, adh_table)
    for i in indexes:
        i = tuple(i)
        # The second keys and locks are switched on purpose so we can compare key1 - key2 and lock1 - lock2
        weighted_diff[i] = max_diff - cell_contact_energy(jk_matrix[0, i[1], i[2]],
                                                          jl_matrix[0, i[1], i[2]],
                                                          jl_matrix[i],
                                                          jk_matrix[i],
                                                          adh_table)
    weighted_count = np.sum(weighted_diff, axis=0)

    return datadf, tau_diff, anc_diffs, anc_diff_counts, weighted_diff, weighted_count


def bit_evolution(datadf: pd.DataFrame):
    bitdf = datadf.drop(columns="tau")
    jkeys = pivot_matrix(bitdf, "jkey_dec")
    jlocks = pivot_matrix(bitdf, "jlock_dec")
    return jkeys, jlocks


def to_bitarray(jseries: pd.Series):
    # Since we are oly computing the hamming distance we dont need to know
    # the actual sizes of the original bitstrings
    length = int(np.log2(jseries.max())) + 1
    return np.array([dec_to_bitarray(x, length) for x in jseries])


def tau_evolution(datadf: pd.DataFrame):
    # Only care about tau for now
    mutdf = datadf.drop(columns=["jkey_dec", "jlock_dec"])
    return pivot_matrix(mutdf, "tau")


def pivot_matrix(df: pd.DataFrame, values):
    pivotdf = pd.pivot(df, index=["id", "foodp"], columns="chem", values=values)
    shape = (df["id"].nunique(), df["foodp"].nunique(), df["chem"].nunique())
    return np.flip(np.reshape(pivotdf.values, shape), 1)


def make_datadf(cells: pd.Series, replicas, mut_rate, mut_std, *sweepargs):
    with tempfile.TemporaryDirectory() as tempdir:
        dirpath = Path(tempdir)
        genomefile = dirpath / "genome.csv"
        write_genomes(cells.to_frame().T, genomefile)
        mutfile = dirpath / "mut.csv"
        mutation_wrapper(genomefile, mutfile, replicas, mut_rate, mut_std)
        sweepfile = dirpath / "sweep.csv"
        sweep_genome_wrapper(mutfile, sweepfile, *sweepargs)
        return pd.read_csv(sweepfile)


# Vectorized hamming distance for bitstrings given as base 10 integers
# This is chatgpt code that i tested and seems to work
# def hamming_distance(b1, b2):
#     # Perform bitwise XOR operation
#     xor_result = b1 ^ b2
#
#     # Count the number of set bits (1s)
#     hamming_dist = 0
#     while xor_result:
#         hamming_dist += xor_result & 1
#         xor_result >>= 1
#
#     return hamming_dist


if __name__ == "__main__":
    main()
