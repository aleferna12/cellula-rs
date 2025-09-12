/*
TODO!:
    - finish IO
        - backup (TEST)
 */
use std::path::PathBuf;
use clap::Parser;
use model::io::io_manager::CONFIG_COPY_PATH;
use model::io::parameters::{Cli, Parameters};
use model::io::parameters::Commands::{Resume, Run};
use model::model::Model;

fn main() -> anyhow::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let cli = Cli::parse();
    
    let parameters = match cli.command {
        Run { config } => Parameters::parse(config),
        Resume { directory, config } => match config {
            Some(config_) => Parameters::parse(config_),
            None => Parameters::parse(PathBuf::from(directory).join(CONFIG_COPY_PATH))
        }
    }?;
    parameters.check_conflicts();
    
    let mut model = Model::initialise_from_parameters(parameters, None)?;
    model.run();
    Ok(())
}
