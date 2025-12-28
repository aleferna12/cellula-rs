//! Contains logic related to the simulation parameters.

use cellulars_lib::constants::CellIndex;
use clap::{Parser, Subcommand};
use serde::{Deserialize, Serialize};
use std::path::Path;
use strum_macros::EnumIter;

static RUN_NOTES: &str = "\
    Model parameters are loaded from a TOML file specified by CONFIG.\n\
    You can also override any parameter from the CONFIG file with environmental variables \
    (use `CPM` as a prefix and `__` as a separator for the parameter section, e.g. `CPM__GENERAL__TIME_STEPS=100`).\n\
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
        /// Path to a TOML file with parameters
        config: String
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
    }
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
    /// Parses parameters from a config file at `path` + env. variables.
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

    /// Checks for conflicting parameters choices and panics if any are found.
    pub fn check_conflicts(&self) -> anyhow::Result<()> {
        #[cfg(not(feature = "fixed-boundary"))]
        if self.pond.enclose && self.pond.neigh_r > 1 {
            anyhow::bail!(
                "`enclose` can only be used with `neigh-r=1`. \
                 If you need an enclosed pond with larger neighbourhoods, enable the `fixed_boundary` feature."
            );
        }
        #[cfg(feature = "fixed-boundary")]
        if !self.pond.enclose {
            anyhow::bail!("`enclose` must be `true` when the `fixed_boundary` feature is enabled")
        }
        Ok(())
    }
}

/// General simulation parameters.
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct GeneralParameters {
    pub time_steps: u32,
    pub seed: Option<u64>
}

/// Parameters determining how a pond is created (see [pond](crate::pond));
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct PondParameters {
    pub width: usize,
    pub height: usize,
    #[serde(default = "param_defaults::false_flag")]
    pub enclose: bool,
    pub neigh_r: u8,
}

/// Parameters specifying how cells are created and behave (see [cell](crate::cell)).
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
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
}

/// Parameters for the cellular automata update algorithm (see [potts](crate::potts)).
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct PottsParameters {
    pub boltz_t: f32,
    pub size_lambda: f32,
    pub chemotaxis_mu: f32,
    pub adhesion: AdhesionParameters
}

/// Parameters used in cell adhesion (see [cellulars_lib::static_adhesion]).
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct AdhesionParameters {
    pub cell_energy: f32,
    pub medium_energy: f32,
    pub solid_energy: f32,
}

/// Parameters used to control IO operations (see [io_manager](crate::io::io_manager)).
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
    pub plot: PlotParameters,
    pub movie: Option<MovieParameters>,
}

/// Parameters used to determine how and when to save data (see [io_manager](crate::io::io_manager)).
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct DataParameters {
    pub cells_period: u32,
    pub lattice_period: u32
}

/// Parameters used to display the real-time movie of the simulation (see [movie_maker](crate::io::movie_maker)).
///
/// Omitting these from the configuration file disables the movie window (same as setting `show` = False).
/// The `movie` feature flag must be on for the movie to be displayed.
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
pub struct MovieParameters {
    #[serde(default = "param_defaults::false_flag")]
    pub show: bool,
    pub width: u32,
    pub height: u32,
    pub frame_period: u32
}

/// Parameters using for plotting (see [plot](crate::io::plot)).
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
    pub migrating_color: String,
    pub dividing_color: String
}


/// A type of plot.
#[derive(Serialize, Deserialize, Clone, EnumIter, Debug)]
#[serde(rename_all = "kebab-case")]
pub enum PlotType {
    /// Cell spin.
    Spin,
    /// Cell center.
    Center,
    /// Cell perceived chemical center.
    ChemCenter,
    /// Cell border.
    Border,
    /// Cell type.
    CellType,
    /// Cell area.
    Area,
    /// Background chemical.
    Chem
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
