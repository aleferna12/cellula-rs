use clap::{Parser, Subcommand};

static RUN_NOTES: &str = "\
    If `LAYOUT` is provided, these specifications are to be followed:\n\t\
        1. Only grayscale images are allowed.\n\t\
        2. Black pixels represent solid objects.\n\t\
        3. White pixels represent the medium.\n\t\
        4. Pixels of any other color are aggregated into square cells, \
           with their area given by the `cell.starting_area` simulation parameter.
    \n\
    Additionally, the image will be resized to match the width and height given in the simulation parameters.
    \n\
    If `TEMPLATE` is provided but `LAYOUT` is not provided, \
    an approximatedly equal number of cells will be initialized at random positions using each template.\
    If both `LAYOUT` and `TEMPLATES` are provided, each color in `LAYOUT` is ordered by their value and assigned \
    to the template at the corresponding index in `TEMPLATES` \
    (there must be at least as many templates as there are colors in the layout).
    \n\
    Model parameters are loaded from a TOML file specified by `CONFIG`.\
    You can also override any parameter from the CONFIG file with environmental variables \
    (use `CPM` as a prefix and `__` as a separator for the parameter section, e.g. `CPM__GENERAL__TIME_STEPS=100`).\
    Use commas to pass parameters that expect lists (e.g. `CPM__IO__PLOT__ORDER=spin,center`).
    \n\
    Documentation for parameters can be found in `model/examples/1_cell.toml`.\
";

/// CLI tool that executes [`Commands`].
#[derive(Parser)]
#[command(version, about)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

/// Commands available to the [`Cli`].
#[derive(Subcommand)]
pub enum Commands {
    /// Start a new run
    #[command(after_long_help = RUN_NOTES)]
    Run {
        /// Path to a TOML file with the simulation parameters
        config: String,
        /// Path to a grayscale PNG file containing the layout of cells to be initialized
        /// (if omitted, cells will be initialized at a random positions)
        #[arg(short, long)]
        layout: Option<String>,
        /// Path to PARQUET file containing cell templates used to initialize cells in the simulation
        /// (if omitted, cells are initialized using the simulation parameters)
        #[arg(short, long)]
        templates: Option<String>

    },
    /// Resume a previous run
    Resume {
        /// Path to the directory of the simulation to be resumed
        directory: String,
        /// Time step from which to restore the data from (if omitted, the last time-step will be used)
        #[arg(short, long)]
        time_step: Option<u32>,
        #[arg(short, long)]
        /// Path to a TOML file with parameters (if omitted, will read parameters from the run's `config.toml` file)
        config: Option<String>
    },

}