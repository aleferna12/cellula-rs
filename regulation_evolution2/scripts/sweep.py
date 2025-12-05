import tempfile
import logging
import subprocess
import pandas as pd
from pathlib import Path
from scripts.data_processing import get_parameter_range, filter_kde, make_density_matrix
from scripts.fileio import write_genomes, parse_cell_data
from scripts.config import PROJECT_DIR

BIN_DIR = PROJECT_DIR / "bin"

logger = logging.getLogger(__name__)


def sweep_genome_wrapper(inputfile,
                         outputfile,
                         min_chem,
                         max_chem,
                         step_chem,
                         min_foodp,
                         max_foodp,
                         step_foodp,
                         mcss,
                         reset):
    args = [BIN_DIR / "sweep_genomes", inputfile, outputfile, min_chem, max_chem, step_chem,
            min_foodp, max_foodp, step_foodp, mcss, str(reset).lower()]
    return subprocess.run([str(arg) for arg in args], capture_output=True, check=True)


def sweep_cells(celldf, parameter_range=None, **sweepkwargs):
    if parameter_range is not None:
        sweepkwargs["min_chem"] = parameter_range.loc["min", "grad_conc"]
        sweepkwargs["max_chem"] = parameter_range.loc["max", "grad_conc"]
        sweepkwargs["min_foodp"] = parameter_range.loc["min", "food"] / 204  # Assumes grn_update_period etc!!
        sweepkwargs["max_foodp"] = parameter_range.loc["max", "food"] / 204

    with tempfile.TemporaryDirectory() as tempdir:
        dirpath = Path(tempdir)
        genomefile = dirpath / "genome.csv"
        write_genomes(celldf, genomefile)
        sweepfile = dirpath / "sweep.csv"
        sweep_genome_wrapper(
            genomefile,
            sweepfile,
            **sweepkwargs
        )
        sweepdf = pd.read_csv(sweepfile)
        # Convenience
        sweepdf["bitstring"] = sweepdf["jkey_dec"].astype(str) + "-" + sweepdf["jlock_dec"].astype(str)
        return sweepdf.set_index(celldf.index[sweepdf["id"]])


def sweep_cell(cellss: pd.Series, parameter_range=None, **sweepkwargs):
    return sweep_cells(cellss.to_frame().T, parameter_range=parameter_range, **sweepkwargs)


def mutation_wrapper(inputfile, outputfile, replicas, mut_rate, mut_std):
    args = [BIN_DIR / "mutation", inputfile, outputfile, replicas, 1, mut_rate, mut_std]
    return subprocess.run([str(arg) for arg in args], capture_output=True)


# This arguments cannot be estimated from range_data
def add_inestimable_sweep_args(parser):
    parser.add_argument("--include-thres",
                        help="How much of the real data should fit within the parameter ranges, expressed "
                             "in the range [0, 1], (default: %(default)s)",
                        default=0.95,
                        type=float)
    parser.add_argument("--step-chem",
                       help="Step chem used for sweeping the genome (default: %(default)s)",
                       default=5,
                       type=float)
    parser.add_argument("--step-foodp",
                       help="Step food parcels (not the same as food!) used for sweeping the genome "
                            "(default: %(default)s)",
                       default=1,
                       type=float)
    parser.add_argument("--mcss",
                       help="See docs for sweep_genome (default: %(default)s)",
                       default=10,
                       type=int)
    parser.add_argument("--reset",
                       help="See docs for sweep_genome (default: %(default)s)",
                       action="store_true")


def add_estimable_sweep_args(parser):
    parser.add_argument("--range-data",
                        help="Data folder with cell data from which the parameter ranges for the "
                             "parameter sweep can be automatically estimated")
    parser.add_argument("--min-chem",
                        help="Min chem used for sweeping the genome (required if 'range-data is not set')",
                        type=float)
    parser.add_argument("--max-chem",
                        help="Max chem used for sweeping the genome (required if 'range-data is not set')",
                        type=float)
    parser.add_argument("--min-foodp",
                        help="Min food parcels (not the same as food!) used for sweeping the genome "
                             "(required if 'range-data is not set')",
                        type=float)
    parser.add_argument("--max-foodp",
                        help="Max food parcels (not the same as food!) used for sweeping the genome "
                             "(required if 'range-data is not set')",
                        type=float)


def add_sweep_args(parser):
    add_estimable_sweep_args(parser)
    add_inestimable_sweep_args(parser)


def parse_sweep_args(args):
    if getattr(args, "range_data", True) is not None:
        sweepkwargs = {
            "step_chem": args.step_chem,
            "step_foodp": args.step_foodp,
            "mcss": args.mcss,
            "reset": args.reset
        }
        return sweepkwargs

    sweepkwargs = {
        "min_chem": args.min_chem,
        "max_chem": args.max_chem,
        "step_chem": args.step_chem,
        "min_foodp": args.min_foodp,
        "max_foodp": args.max_foodp,
        "step_foodp": args.step_foodp,
        "mcss": args.mcss,
        "reset": args.reset
    }
    if None in sweepkwargs.values():
        raise Exception("either 'range-data' or min and max arguments must be set")

    return sweepkwargs
