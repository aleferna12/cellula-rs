use clap::Parser;

// TODO: implement Display
#[derive(Parser, Debug, Default)]
#[command(version, about, long_about = None)]
pub struct Parameters {
    #[arg(long, default_value_t = 1_000)]
    pub time_steps: u32,
    #[arg(long, default_value_t = 100)]
    pub n_cells: u16,
    #[arg(long, default_value_t = 50)]
    pub cell_area: u32,
    #[arg(long, default_value_t = 50)]
    pub cell_target_area: u32,
    #[arg(long, default_value_t = 100)]
    pub width: usize,
    #[arg(long, default_value_t = 100)]
    pub height: usize,
    #[arg(long, default_value_t = 0)]
    pub seed: u64,
    #[arg(long, default_value_t = 1)]
    pub neigh_r: u8,
    #[arg(long, default_value_t = 16.)]
    pub boltz_t: f64,
    #[arg(long, default_value_t = 1.)]
    pub size_lambda: f64,
    #[arg(long, default_value_t = 16.)]
    pub cell_energy: f64,
    #[arg(long, default_value_t = 16.)]
    pub med_energy: f64,
    #[arg(long, default_value_t = 16.)]
    pub solid_energy: f64
}