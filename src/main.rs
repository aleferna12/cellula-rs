use std::error::Error;
use clap::Parser;
use evo_cpm::model::Model;
use evo_cpm::parameters::Parameters;
use evo_cpm::io::welcome;
/*
TODO!: 
    - read params from config file
    - finish IO
        - cell info
        - backup
 */

fn main() -> Result<(), Box<dyn Error>> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    
    let parameters = Parameters::parse();
    parameters.check_conflicts();
    welcome(&parameters);
    
    let mut model = Model::new(parameters);
    model.setup()?;
    model.run(model.parameters.time_steps)?;
    Ok(())
}
