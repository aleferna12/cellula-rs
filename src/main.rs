/*
TODO!: 
    - workout what to keep from WrappedPos 
    - finish IO
        - cell info
        - backup
    - seed is not saved with the model when its 0
    - look into anyhow, I think the backtracing would be very useful
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
    
    let mut model = Model::try_from(parameters)?;
    model.run(model.parameters().general.time_steps);
    Ok(())
}
