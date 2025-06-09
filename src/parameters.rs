use clap::Parser;

#[derive(Parser, Debug)]
#[command(name = "EvoCPM")]
#[command(version, about, long_about = None)]
pub struct Parameters {
    #[arg(long, default_value_t = 100)]
    pub width: usize,
    #[arg(long, default_value_t = 100)]
    pub height: usize,
}