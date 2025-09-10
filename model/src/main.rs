/*
TODO!:
    - finish IO
        - backup (TEST)
 */
use clap::Parser;
use model::io::parameters::{Cli, Parameters};
use model::model::Model;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let cli = Cli::parse();
    
    let parameters = Parameters::parse(&cli.config)?;
    parameters.check_conflicts();
    
    let mut model = Model::initialise_from_parameters(parameters)?;
    model.run();
    Ok(())
}
