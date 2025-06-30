use std::path::Path;
use clap::Parser;
use serde::{Deserialize, Serialize};

static CLI_NOTES: &str = "\
    Model parameters are loaded from a TOML file specified by CONFIG.\n\
    You can also override any parameter from the CONFIG file with environmental variables \
    (use `__` for the parameter section, e.g. `GENERAL__TIME_STEPS`).\n\
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
    pub io: IoParameters,
    pub movie: MovieParameters
}

impl Parameters {
    pub fn from_file(path: impl AsRef<Path>) -> Result<Parameters, config::ConfigError> {
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
            ).build()?
                .try_deserialize()?;
        Ok(params)
    }
    
    pub fn check_conflicts(&self) {
        if self.environment.enclose && self.environment.neigh_r > 1 {
            log::warn!("`enclose` can only be used when `neigh-r` == 1 by default");
            log::warn!("You can circumvent this issue by changing the `Boundary` type in `Environment` \
                        from `UnsafePeriodicBoundary` to `FixedBoundary`");
        }
    }
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub struct GeneralParameters {
    pub time_steps: u32,
    #[serde(default = "param_defaults::seed")]
    pub seed: u64,
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub struct EnvironmentParameters {
    pub width: usize,
    pub height: usize,
    #[serde(default = "param_defaults::enclose")]
    pub enclose: bool,
    pub n_cells: u32,
    pub cell_start_area: u32,
    pub neigh_r: u8,
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub struct CellularAutomataParameters {
    pub cell_target_area: u32,
    pub boltz_t: f32,
    pub size_lambda: f32,
    pub adhesion: AdhesionParameters
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(tag = "type", rename_all = "kebab-case")]
pub enum AdhesionParameters {
    #[serde(rename_all = "kebab-case")]
    StaticAdhesion {
        cell_energy: f32,
        medium_energy: f32,
        solid_energy: f32,
    }
}

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub struct IoParameters {
    pub outdir: String,
    #[serde(default = "param_defaults::replace_outdir")]
    pub replace_outdir: bool,
    pub image_period: u32,
    #[serde(default = "param_defaults::image_format")]
    pub image_format: String,
}

// TODO! docs
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub struct MovieParameters {
    #[serde(default = "param_defaults::show")]
    pub show: bool,
    pub width: u32,
    pub height: u32,
    pub frame_period: u32
}

// This is a workaround while https://github.com/serde-rs/serde/issues/368 is pending
mod param_defaults {
    pub fn seed() -> u64 { 0 }
    pub fn replace_outdir() -> bool { false }
    pub fn image_format() -> String { "webp".to_string() }
    pub fn enclose() -> bool { false }
    pub fn show() -> bool { false }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_from_file() -> Result<(), config::ConfigError> {
        Parameters::from_file("examples/64_cells.toml")?;
        Ok(())
    }
}
