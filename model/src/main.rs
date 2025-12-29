//! Entry point for the main binary.
//!
//! Check [Cli] for the available commands (or run with `--help`).

/*
TODO!:
    - profile data write, something is taking awfully long (prob lattice write)
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
        Run {
            config,
            layout: maybe_layout,
            templates: maybe_templates
        } => {
            let params = Parameters::parse(config)?;
            // TODO!: implement templates
            match maybe_layout {
                None => Model::new_from_parameters(params)?,
                Some(layout) => Model::new_from_layout(params, layout)?
            }
        },
        Resume {
            directory,
            config: maybe_config,
            time_step: maybe_time_step
        } => {
            let params = match maybe_config {
                Some(config) => Parameters::parse(config),
                None => Parameters::parse(IoManager::resolve_parameters_path(&directory))
            }?;
            let time_step = maybe_time_step.unwrap_or(IoManager::find_last_time_step(&directory)?);
            Model::new_from_backup(params, directory, time_step)?
        }
    };
    model.run();
    model.goodbye();
    Ok(())
}