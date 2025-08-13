use std::error::Error;
use rand::distr::{Distribution, Uniform};
use rand::{RngCore, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use crate::adhesion::{ClonalAdhesion, StaticAdhesion};
use crate::cell::Cell;
use crate::cell_container::CellContainer;
use crate::cellular_automata::CellularAutomata;
use crate::constants::{BoundaryType, NeighbourhoodType};
use crate::environment::{Environment, LatticeEntity};
use crate::genome::Grn;
use crate::io::io_manager::IoManager;
use crate::io::movie_maker::MovieMaker;
use crate::io::parameters::Parameters;
use crate::pond::Pond;
use crate::positional::rect::Rect;
use crate::space::Space;

pub struct Model {
    pub ponds: Vec<Pond>,
    pub io: IoManager,
    time_steps: u32
}

impl Model {
    pub fn initialise_from_parameters(mut parameters: Parameters) -> Result<Model, Box<dyn Error>> {
        log::info!("Initialising model");
        // TOML doesnt support large u64s so we use a u32 seed
        let seed = parameters.general.seed.unwrap_or(Xoshiro256StarStar::from_os_rng().next_u32() as u64);
        parameters.general.seed = seed.into();
        let mut rng = Xoshiro256StarStar::seed_from_u64(seed);

        let io = IoManager::new(
            &parameters.io.outdir,
            parameters.io.image_period,
            parameters.io.image_format.clone(),
            parameters.io.plot.clone(),
            if parameters.io.movie.show {
                match MovieMaker::new(
                    parameters.io.movie.width,
                    parameters.io.movie.height,
                    parameters.io.movie.frame_period
                ) {
                    Ok(mm) => {
                        log::info!("Creating window for real-time movie display");
                        Some(mm)
                    },
                    Err(e) => {
                        log::warn!("Failed to initialise movie maker with error `{e}`");
                        None
                    }
                }
            } else {
                None
            }
        );

        log::info!("Creating output directories and copy of parameter file");
        io.create_directories(parameters.io.replace_outdir)?;
        io.create_parameters_file(&parameters)?;

        let mut ponds = vec![];
        // TODO: if we make everything clonable then that helps here and also in model_bench
        //  (takes less time to reinitialise everything which means more samples)
        for pond_i in 0..parameters.general.n_environments {
            log::info!("Making pond #{pond_i}");
            let mut env = Environment::new(
                parameters.environment.update_period,
                parameters.environment.cell_search_radius,
                parameters.environment.max_cells,
                parameters.environment.enclose,
                CellContainer::new(
                    parameters.cell.target_area,
                    parameters.cell.divide,
                    parameters.cell.migrate,
                ),
                Space::new(
                    BoundaryType::new(Rect::new(
                        (0., 0.).into(),
                        (parameters.environment.width as f32, parameters.environment.height as f32).into(),
                    ))
                )?,
                NeighbourhoodType::new(parameters.environment.neigh_r)
            );

            log::info!("Creating cells");
            let cell_side = (parameters.cell.starting_area as f32).sqrt() as usize;
            for _ in 0..parameters.environment.starting_cells {
                let pos = env.space.cell_lattice.random_pos(&mut rng);
                env.spawn_rect_cell(
                    Rect::new(
                        pos,
                        (pos.x + cell_side, pos.y + cell_side).into()
                    ),
                    Cell::new_empty(
                        parameters.cell.target_area,
                        parameters.cell.div_area,
                        Grn::new(
                            [1. / env.height() as f32],
                            parameters.cell.n_regulatory_genes,
                            parameters.cell.mutation_rate,
                            parameters.cell.mutation_std,
                            || Uniform::new(-1., 1.).unwrap().sample(&mut rng)
                        )
                    )
                );
            }
            log::info!(
                "Created {} out of the {} cells requested",
                env.cells.n_cells(),
                parameters.environment.starting_cells
            );

            let ca= CellularAutomata::new(
                parameters.cellular_automata.boltz_t,
                parameters.cellular_automata.size_lambda,
                parameters.cellular_automata.chemotaxis_mu,
                ClonalAdhesion::new(
                    parameters.environment.max_cells + LatticeEntity::first_cell_spin(),
                    StaticAdhesion {
                        cell_energy: parameters.cellular_automata.adhesion.cell_energy,
                        medium_energy: parameters.cellular_automata.adhesion.medium_energy,
                        solid_energy: parameters.cellular_automata.adhesion.solid_energy,
                    }
                )
            );
            ponds.push(Pond::new(env, ca, rng.clone()));
        }

        Ok(Self { ponds, io, time_steps: parameters.general.time_steps })
    }

    pub fn run_for(&mut self, time_steps: u32) {
        for time_step in 0..=time_steps {
            let saved = self.io.image_io(
                time_step,
                &self.ponds
            );
            if let Err(e) = saved {
                log::warn!("Failed to save image at time step {time_step} with error `{e}`")
            }
            for pond in &mut self.ponds {
                pond.step();
            }
        }
    }

    pub fn run(&mut self) {
        self.run_for(self.time_steps);
    }
}

#[cfg(test)]
mod tests {
    use rand::{Rng, SeedableRng};
    use rand_xoshiro::Xoshiro256StarStar;

    #[test]
    fn test_seed() {
        let mut rng = Xoshiro256StarStar::seed_from_u64(1241254152);
        let s = (0..50)
            .map(|_| rng.random_range(0..9).to_string())
            .collect::<Vec<_>>()
            .join("");
        let res = "15515320360704325727185856564110164830043067488704";
        assert_eq!(res, s);
    }
}