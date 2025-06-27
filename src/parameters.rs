use clap::Parser;
use serde::{Deserialize, Serialize};

static CLI_NOTES: &str = "\
    Model parameters are loaded from a TOML file specified by CONFIG.\n\
    You can also override any parameter from the CONFIG file with environmental variables \
    (use the `CPM` prefix and UPPER_SNAKE_CASE, e.g. `CPM_TIME_STEPS`).\n\
    Documentation for parameters can be found in `examples/64_cells.toml`.\n\
";

#[derive(Parser)]
#[command(version, about, after_long_help = CLI_NOTES)]
pub struct Cli {
    #[arg(help = "Path to TOML file storing the model parameters")]
    pub config: String
}

/// Parameters for the model.
///
/// Documentation for each parameter is in `examples/64_cells.toml`
#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "kebab-case")]
pub struct Parameters {
    pub time_steps: u32,
    pub seed: u64,
    pub outdir: String,
    pub replace_outdir: bool,
    pub image_period: u32,

    pub width: usize,
    pub height: usize,
    pub enclose: bool,
    pub n_cells: u32,

    pub cell_start_area: u32,
    pub cell_target_area: u32,
    pub boltz_t: f32,
    pub neigh_r: u8,
    pub size_lambda: f32,
    pub cell_energy: f32,
    pub medium_energy: f32,
    pub solid_energy: f32,
}

impl Parameters {
    pub fn check_conflicts(&self) {
        if self.enclose && self.neigh_r > 1 {
            log::warn!("`--enclose` can only be used when `--neigh-r` == 1 by default");
            log::warn!("You can circumvent this issue by changing the `Boundary` type in `Environment` \
                   from `UnsafePeriodicBoundary` to `FixedBoundary`");
        }
    }
}
