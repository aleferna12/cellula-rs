/*
TODO!:
    - finish IO
        - cell info
        - backup
    - add builder structs for the long constructors (Env, IoManager, CellCont...)
    - generalise remaining modules:
        - Space
        - CA
        - Pond (at that point also move functions that modify environment out of environment)
        - At the end, CA should be responsible for executing the steps (only),
            while Pond calls the steps and manages creation/killing of cells
 */
use clap::Parser;
use model::io::parameters::{Cli, Parameters};
use std::error::Error;
use model::model::Model;

fn main() -> Result<(), Box<dyn Error>> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let cli = Cli::parse();
    
    let parameters = Parameters::parse(&cli.config)?;
    parameters.check_conflicts();
    
    let mut model = Model::initialise_from_parameters(parameters)?;
    model.run();
    Ok(())
}
