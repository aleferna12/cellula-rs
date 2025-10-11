use cellulars_lib::constants::CellIndex;
use clap::{Parser, Subcommand};
use serde::{Deserialize, Serialize};
use std::path::Path;
use strum_macros::EnumIter;

static RUN_NOTES: &str = "\
    Model parameters are loaded from a TOML file specified by CONFIG.\n\
    You can also override any parameter from the CONFIG file with environmental variables \
    (use `__` for the parameter section, e.g. `GENERAL__TIME_STEPS`).\n\
    Use commas to pass parameters that expect lists (e.g. `IO__PLOT__ORDER=spin,center`).
    \n\
    Documentation for parameters can be found in `examples/64_cells.toml`.\n\
";

static RESUME_NOTES: &str = "\
    Model parameters can be specified via CONFIG, or by setting environmental variables \
    (see help of the `run` subcommand).\n\
    If CONFIG is not specified and no environmental variables are set, \
    the simulation runs with its original parameters.
";

#[derive(Parser)]
#[command(version, about)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Start a new run
    #[command(after_long_help = RUN_NOTES)]
    Run {
        /// Path to a TOML file with parameters
        config: String
    },
    /// Resume a previous run
    #[command(after_help = RESUME_NOTES)]
    Resume {
        /// Path to the directory of the simulation to be resumed
        directory: String,
        // TODO!: make optional and find last time_step
        /// Time step from which to restore the data from
        time_step: u32,
        /// Path to a TOML file with parameters
        config: Option<String>
    }
}

// When you add parameters, dont forget to document them (and their defaults)
/// Parameters for the model.
///
/// Documentation for each parameter is in `examples/64_cells.toml`
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub struct Parameters {
    pub general: GeneralParameters,
    pub pond: PondParameters,
    pub cell: CellParameters,
    pub potts: PottsParameters,
    pub io: IoParameters
}

impl Parameters {
    pub fn parse(path: impl AsRef<Path>) -> Result<Parameters, config::ConfigError> {
        let path = path.as_ref();
        log::info!("Reading parameters from {} and environmental variables", path.display());
        let params = config::Config::builder()
            .add_source(
                config::File::from(path)
            ).add_source(
                // Converts an env CPM_TIME_STEPS to time-steps
                config::Environment::default()
                    .separator("__")
                    .convert_case(config::Case::Kebab)
                    .list_separator(",")
                    .with_list_parse_key("io.plot.order")
                    .try_parsing(true)
            ).build()?.try_deserialize()?;
        Ok(params)
    }
    
    pub fn check_conflicts(&self) {
        if self.pond.enclose && self.pond.neigh_r > 1 {
            log::warn!("`enclose` can only be used when `neigh-r=1` by default");
            log::warn!("You can circumvent this issue by changing `LatticeBoundaryType` in `Model` \
                        from `UnsafePeriodicBoundary` to `FixedBoundary`");
        }
    }
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub struct GeneralParameters {
    pub time_steps: u32,
    pub seed: Option<u64>,
    pub dispersion_period: u32
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub struct PondParameters {
    pub width: usize,
    pub height: usize,
    #[serde(default = "param_defaults::false_flag")]
    pub enclose: bool,
    pub neigh_r: u8,
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub struct CellParameters {
    pub starting_cells: CellIndex,
    pub max_cells: CellIndex,
    pub search_radius: f32,
    pub starting_area: u32,
    pub target_area: u32,
    pub div_area: u32,
    #[serde(default = "param_defaults::true_flag")]
    pub divide: bool,
    #[serde(default = "param_defaults::true_flag")]
    pub migrate: bool,
    pub update_period: u32,
    pub genome: GenomeParameters,
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub struct GenomeParameters {
    pub n_regulatory: usize,
    pub mutation_rate: f32,
    pub mutation_std: f32,
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub struct PottsParameters {
    pub boltz_t: f32,
    pub size_lambda: f32,
    pub chemotaxis_mu: f32,
    pub act_max: u32,
    pub act_lambda: f32,
    pub adhesion: ClonalAdhesionParameters
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub struct ClonalAdhesionParameters {
    pub clone_energy: f32,
    pub cell_energy: f32,
    pub medium_energy: f32,
    pub solid_energy: f32,
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
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
#[serde(rename_all = "kebab-case")]
pub struct DataParameters {
    pub cells_period: u32,
    pub genomes_period: u32,
    pub clones_period: u32,
    pub lattice_period: u32
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub struct MovieParameters {
    #[serde(default = "param_defaults::false_flag")]
    pub show: bool,
    pub width: u32,
    pub height: u32,
    pub frame_period: u32
}

// We flatten the parameters here to allow order to be an env variable
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub struct PlotParameters {
    pub order: Vec<PlotType>, 
    pub solid_color: String,
    pub medium_color: Option<String>,
    pub center_color: String,
    pub chem_center_color: String,
    pub clones_color: String,
    #[serde(default = "param_defaults::false_flag")]
    pub all_clones: bool,
    pub border_color: String,
    pub area_min_color: String,
    pub area_max_color: String,
    pub chem_min_color: String,
    pub chem_max_color: String,
    pub act_min_color: String,
    pub act_max_color: String,
    pub migrating_color: String,
    pub dividing_color: String
}

#[derive(Serialize, Deserialize, Clone, EnumIter, Debug)]
#[serde(rename_all = "kebab-case")]
pub enum PlotType {
    Spin,
    Center,
    ChemCenter,
    Clones,
    Border,
    CellType,
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
    fn test_parse() -> Result<(), config::ConfigError> {
        Parameters::parse("examples/64_cells.toml")?;
        Ok(())
    }
}
