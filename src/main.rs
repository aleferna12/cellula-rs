/*
TODO!:
    - finish IO
        - cell info
        - backup
    - add builder structs for the long constructors (Env, IoManager, CellCont...)
 */
use clap::Parser;
use evo_cpm::io::parameters::{Cli, Parameters};
use std::error::Error;
use evo_cpm::model::Model;

fn main() -> Result<(), Box<dyn Error>> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let cli = Cli::parse();
    
    let parameters = Parameters::parse(&cli.config)?;
    parameters.check_conflicts();
    
    let mut model = Model::initialise_from_parameters(parameters)?;
    model.run();
    Ok(())
}
