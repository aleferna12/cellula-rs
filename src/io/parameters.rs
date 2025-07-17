use crate::constants::Spin;
use clap::Parser;
use serde::{Deserialize, Serialize};
use std::path::Path;

static CLI_NOTES: &str = "\
    Model parameters are loaded from a TOML file specified by CONFIG.\n\
    You can also override any parameter from the CONFIG file with environmental variables \
    (use `__` for the parameter section, e.g. `GENERAL__TIME_STEPS`).\n\
    Use commas to pass parameters that expect lists (e.g. `IO__PLOT__ORDER=spin,center`).
    \n\
    Documentation for parameters can be found in `examples/64_cells.toml`.\n\
";

#[derive(Parser)]
#[command(version, about, after_long_help = CLI_NOTES)]
pub struct Cli {
    #[arg(help = "Path to TOML file storing the model parameters")]
    pub config: String
}

// When you add parameters, dont forgot to document them (and their defaults)
/// Parameters for the model.
///
/// Documentation for each parameter is in `examples/64_cells.toml`
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub struct Parameters {
    pub general: GeneralParameters,
    pub environment: EnvironmentParameters,
    pub cellular_automata: CellularAutomataParameters,
    pub io: IoParameters
}

impl Parameters {
    pub fn parse(path: impl AsRef<Path>) -> Result<Parameters, config::ConfigError> {
        let path = path.as_ref();
        log::info!("Reading parameters from `{}` and environment", path.display());
        let params = config::Config::builder()
            .add_source(
                config::File::from(path)
            ).add_source(
                // Converts an env CPM_TIME_STEPS to time-steps
                config::Environment::default()
                    .separator("__")
                    .convert_case(config::Case::Kebab)
                    .list_separator(",")
                    .with_list_parse_key("io.plots.order")
                    .try_parsing(true)
            ).build()?
                .try_deserialize()?;
        Ok(params)
    }
    
    pub fn check_conflicts(&self) {
        if self.environment.enclose && self.environment.neigh_r > 1 {
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
    #[serde(default = "param_defaults::zero")]
    pub seed: u64,
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub struct EnvironmentParameters {
    pub width: usize,
    pub height: usize,
    #[serde(default = "param_defaults::false_flag")]
    pub enclose: bool,
    pub neigh_r: u8,
    pub starting_cells: Spin,
    pub max_cells: Spin,
    pub cell_start_area: u32,
    pub cell_search_radius: f32,
    pub update_period: u32,
    pub cell: CellParameters
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub struct CellParameters {
    pub target_area: u32,
    pub div_area: u32,
    #[serde(default = "param_defaults::true_flag")]
    pub divide: bool,
    // TODO: change to true when migration is properly implemented (its currently an unstable feature)
    #[serde(default = "param_defaults::false_flag")]
    pub migrate: bool
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub struct CellularAutomataParameters {
    pub boltz_t: f32,
    pub size_lambda: f32,
    pub chemotaxis_mu: f32,
    pub adhesion: StaticAdhesionParameters
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub struct StaticAdhesionParameters {
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
    pub movie: MovieParameters,
    pub plots: PlotParameters
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
    pub clones_color: String,
    #[serde(default = "param_defaults::false_flag")]
    pub all_clones: bool,
    pub border_color: String,
    pub area_min_color: String,
    pub area_max_color: String,
    pub light_min_color: String,
    pub light_max_color: String,
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub enum PlotType {
    Spin,
    Center,
    Clones,
    Area,
    Border,
    Light
}

// This is a workaround while https://github.com/serde-rs/serde/issues/368 is pending
mod param_defaults {
    pub fn zero() -> u64 { 0 }
    pub fn false_flag() -> bool { false }
    pub fn true_flag() -> bool { true }
    pub fn webp() -> String { "webp".to_string() }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_file() -> Result<(), config::ConfigError> {
        Parameters::parse("examples/64_cells.toml")?;
        Ok(())
    }
}
