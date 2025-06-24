use clap::Parser;
use log;

// TODO: implement Display
#[derive(Parser, Debug, Default)]
#[command(version, about, long_about = None)]
pub struct Parameters {
    // General parameters
    /// How many MCS the simulation will run for.
    #[arg(long, default_value_t = 1_001)]
    pub time_steps: u32,
    /// Seed for the RNG.
    ///
    /// Use `0` to pick a random seed.
    #[arg(long, default_value_t = 0)]
    pub seed: u64,
    
    // IO parameters
    /// Output directory for files.
    #[arg(long, default_value = "./out")]
    pub outdir: String,
    /// Flag used to delete contents of `outdir` if it exists.
    #[arg(long, action)]
    pub replace_outdir: bool,
    /// Period (in MCS) with which to save an image of the simulation.
    #[arg(long, default_value_t = 1_000)]
    pub image_period: u32,
    
    // Environment parameters
    /// Width of the lattice.
    #[arg(long, default_value_t = 100)]
    pub width: usize,
    /// Height of the lattice.
    #[arg(long, default_value_t = 100)]
    pub height: usize,
    /// Whether to create a solid border around the cell lattice.
    ///
    /// <div class="warning">Should not be used with `neigh_r` > 1.</div>
    #[arg(long, action)]
    pub enclose: bool,
    /// Number of starting cells in the environment.
    #[arg(long, default_value_t = 100)]
    pub n_cells: u32,
    
    // Cell parameters
    /// Starting cell area.
    #[arg(long, default_value_t = 50)]
    pub cell_start_area: u32,
    /// Target cell area.
    #[arg(long, default_value_t = 50)]
    pub cell_target_area: u32,
    
    // CA parameters
    /// Boltzmann temperature.
    ///
    /// Increases likelihood of unfavourable copy attempts being accepted.
    #[arg(long, default_value_t = 16.)]
    pub boltz_t: f32,
    /// Radius of the neighbourhood.
    #[arg(long, default_value_t = 1)]
    pub neigh_r: u8,
    /// Scaler constant for area deviation penalties.
    #[arg(long, default_value_t = 1.)]
    pub size_lambda: f32,
    /// Adhesion energy at cell-cell interfaces.
    #[arg(long, default_value_t = 16.)]
    pub cell_energy: f32,
    /// Adhesion energy at cell-medium interfaces.
    #[arg(long, default_value_t = 16.)]
    pub medium_energy: f32,
    /// Adhesion energy at cell-solid interfaces.
    #[arg(long, default_value_t = 16.)]
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

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_parse_cli() {
        assert!(Parameters::try_parse_from([""]).is_ok());
    }
}