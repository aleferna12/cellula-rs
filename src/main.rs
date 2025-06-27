/*
TODO!: 
    - finish IO
        - cell info
        - backup
 */
use std::error::Error;
use clap::Parser;
use evo_cpm::model::Model;
use evo_cpm::parameters::{Cli, Parameters};

fn main() -> Result<(), Box<dyn Error>> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let cli = Cli::parse();
    
    let parameters = Parameters::from_file(&cli.config)?;
    parameters.check_conflicts();
    
    let mut model = Model::new(parameters);
    model.setup()?;
    model.run(model.parameters().time_steps);
    Ok(())
}
