/*
TODO!: 
    - profile, there was a 6% performance loss since the last commit 
    - finish IO
        - cell info
        - backup
    - seed is not saved with the model when its 0
    - look into anyhow, I think the backtracing would be very useful
    - add builder structs for the long constructors (Env, IoManager, CellCont...)
 */
use clap::Parser;
use evo_cpm::io::parameters::{Cli, Parameters};
use evo_cpm::model::Model;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let cli = Cli::parse();
    
    let parameters = Parameters::parse(&cli.config)?;
    parameters.check_conflicts();
    
    let time_steps = parameters.general.time_steps;
    let mut model = Model::try_from(parameters)?;
    model.run(time_steps);
    Ok(())
}
