/*
TODO!: 
    - finish IO
        - cell info
        - backup
    - save_image is being called twice when image_period and frame_period match
    - seed is not saved with the model when its 0
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
    
    let mut model = Model::try_from(parameters)?;
    model.run(model.parameters().general.time_steps);
    Ok(())
}
