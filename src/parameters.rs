use clap::Parser;

// TODO: implement Display
#[derive(Parser, Debug, Default)]
#[command(version, about, long_about = None)]
pub struct Parameters {
    #[arg(long, default_value_t = 100_000)]
    pub time_steps: u32,
    #[arg(long, default_value_t = 100)]
    pub width: usize,
    #[arg(long, default_value_t = 100)]
    pub height: usize,
    #[arg(long, default_value_t = 0)]
    pub seed: u64,
    #[arg(long, default_value_t = 1)]
    pub neigh_r: u8,
    #[arg(long, default_value_t = 12f32)]
    pub boltz_t: f32,
    #[arg(long, default_value_t = 1f32)]
    pub size_lambda: f32,
    #[arg(long, default_value_t = 50)]
    pub target_area: u32,
}