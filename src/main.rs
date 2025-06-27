/*
TODO!: 
    - read params from config file
    - finish IO
        - cell info
        - backup
 */
use std::error::Error;
use clap::Parser;
use evo_cpm::model::Model;
use evo_cpm::parameters::Cli;
use evo_cpm::io::read_config;

fn main() -> Result<(), Box<dyn Error>> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let cli = Cli::parse();
    
    let parameters = read_config(&cli.config)?;
    parameters.check_conflicts();
    
    let mut model = Model::new(parameters);
    model.setup()?;
    model.run(model.parameters().time_steps);
    Ok(())
}
