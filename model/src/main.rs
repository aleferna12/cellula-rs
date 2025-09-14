/*
TODO!:
    - finish IO
        - backup (TEST)
    - strip down of non essential features
 */
use clap::Parser;
use model::io::io_manager::IoManager;
use model::io::parameters::Commands::{Resume, Run};
use model::io::parameters::{Cli, Parameters};
use model::model::Model;

fn main() -> anyhow::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let cli = Cli::parse();
    let mut model = match cli.command {
        Run { config } => {
            let params = Parameters::parse(config)?;
            Model::initialise_from_parameters(params)?
        },
        Resume { directory, config, time_step } => {
            let params = match config {
                Some(config_) => Parameters::parse(config_),
                None => Parameters::parse(IoManager::resolve_parameters_path(&directory))
            }?;
            Model::initialise_from_backup(params, directory, time_step)?
        }
    };
    model.run();
    Ok(())
}
