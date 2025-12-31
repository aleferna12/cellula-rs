use cellulars_lib::constants::CellIndex;
use clap::{Parser, Subcommand};
use serde::{Deserialize, Serialize};
use std::path::Path;
use strum_macros::EnumIter;

static RUN_NOTES: &str = "\
    Model parameters are loaded from a TOML file specified by CONFIG.\n\
    You can also override any parameter from the CONFIG file with environmental variables \
    (use `CPM` as a prefix and `__` as a separator for the parameter section, e.g. `CPM__GENERAL__TIME_STEPS`).\n\
    Use commas to pass parameters that expect lists (e.g. `CPM__IO__PLOT__ORDER=spin,center`).
    \n\
    Documentation for parameters can be found in `model/examples/64_cells.toml`.\n\
";

/// CLI tool that executes [Commands].
#[derive(Parser)]
#[command(version, about)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

/// Commands available to the [Cli].
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


// When you add parameters, dont forget to document them (and their defaults)
/// Parameters for the model.
///
/// Documentation for each parameter is in `examples/64_cells.toml`
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct Parameters {
    pub general: GeneralParameters,
    pub pond: PondParameters,
    pub cell: CellParameters,
    pub potts: PottsParameters,
    pub io: IoParameters
}

impl Parameters {
    pub fn parse(path: impl AsRef<Path>) -> anyhow::Result<Parameters> {
        let path = path.as_ref();
        log::info!("Reading parameters from {} and environmental variables", path.display());
        let params: Parameters = config::Config::builder()
            .add_source(
                config::File::from(path)
            ).add_source(
                // Converts an env CPM_TIME_STEPS to time-steps
                config::Environment::default()
                    .prefix("CPM")
                    .prefix_separator("__")
                    .separator("__")
                    .convert_case(config::Case::Kebab)
                    .list_separator(",")
                    .with_list_parse_key("io.plot.order")
                    .try_parsing(true)
            ).build()?
            .try_deserialize()?;
        params.check_conflicts()?;
        Ok(params)
    }
    
    pub fn check_conflicts(&self) -> anyhow::Result<()> {
        #[cfg(not(feature = "fixed-boundary"))]
        if self.pond.enclose && self.pond.neigh_r > 1 {
            anyhow::bail!(
                "`enclose` can only be used with `neigh-r=1`. \
                 If you need an enclosed pond with larger neighbourhoods, enable the `fixed-boundary` feature."
            );
        }
        #[cfg(feature = "fixed-boundary")]
        if !self.pond.enclose {
            anyhow::bail!("`enclose` must be `true` when the `fixed-boundary` feature is enabled")
        }
        #[cfg(not(feature = "static-adhesion"))]
        if self.cell.genome.length < 1 || self.cell.genome.length > 64 {
            anyhow::bail!("`cell.genome.length` must be between 1 and 64")
        }
        Ok(())
    }
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct GeneralParameters {
    pub time_steps: u32,
    pub seed: Option<u64>
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct PondParameters {
    pub width: usize,
    pub height: usize,
    #[serde(default = "param_defaults::false_flag")]
    pub enclose: bool,
    pub neigh_r: u8,
    pub season_duration: u32,
    pub half_fitness: f32,
    pub reproduction_steps: u32
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct CellParameters {
    pub starting_cells: CellIndex,
    pub max_cells: CellIndex,
    pub search_radius: f32,
    pub starting_area: u32,
    pub target_area: u32,
    pub target_perimeter: u32,
    #[serde(default = "param_defaults::true_flag")]
    pub divide: bool,
    pub genome: GenomeParameters,
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub struct GenomeParameters {
    pub mutation_rate: f32,
    pub length: u8,
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct PottsParameters {
    pub boltz_t: f32,
    pub size_lambda: f32,
    pub perimeter_lambda: f32,
    pub chemotaxis_min: f32,
    pub act_max: u32,
    pub act_lambda: f32,
    pub adhesion: AdhesionParameters,
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct AdhesionParameters {
    pub cell_energy: f32,
    pub medium_energy: f32,
    pub solid_energy: f32,
    pub gene_energy: f32
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct IoParameters {
    pub outdir: String,
    #[serde(default = "param_defaults::false_flag")]
    pub replace_outdir: bool,
    pub image_period: u32,
    #[serde(default = "param_defaults::webp")]
    pub image_format: String,
    pub info_period: u32,
    pub data: DataParameters,
    pub movie: MovieParameters,
    pub plot: PlotParameters
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct DataParameters {
    pub cells_period: u32,
    pub lattice_period: u32
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct MovieParameters {
    #[serde(default = "param_defaults::false_flag")]
    pub show: bool,
    pub width: u32,
    pub height: u32,
    pub frame_period: u32
}

// We flatten the parameters here to allow order to be an env variable
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct PlotParameters {
    pub order: Box<[PlotType]>, 
    pub solid_color: String,
    pub medium_color: Option<String>,
    pub center_color: String,
    pub chem_center_color: String,
    pub border_color: String,
    pub area_min_color: String,
    pub area_max_color: String,
    pub chem_min_color: String,
    pub chem_max_color: String,
    pub act_min_color: String,
    pub act_max_color: String
}

#[derive(Serialize, Deserialize, Clone, EnumIter, Debug)]
#[serde(rename_all = "kebab-case")]
pub enum PlotType {
    Spin,
    Center,
    ChemCenter,
    Border,
    Area,
    Chem,
    Act
}

// This is a workaround while https://github.com/serde-rs/serde/issues/368 is pending
mod param_defaults {
    pub fn false_flag() -> bool { false }
    pub fn true_flag() -> bool { true }
    pub fn webp() -> String { "webp".to_string() }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse() -> anyhow::Result<()> {
        Parameters::parse("examples/64_cells.toml")?;
        Ok(())
    }
}
