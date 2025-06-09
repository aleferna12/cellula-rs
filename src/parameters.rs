use std::fmt::{Display, Formatter};
use clap::Parser;

#[derive(Parser, Debug, Default)]
#[command(version, about, long_about = None)]
pub struct Parameters {
    #[arg(long, default_value_t = 100)]
    pub width: usize,
    #[arg(long, default_value_t = 100)]
    pub height: usize,
    #[arg(long, default_value_t = 0)]
    pub seed: u64
}

impl Display for Parameters {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        todo!()
    }
}