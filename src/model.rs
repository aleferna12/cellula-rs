use crate::adhesion::{ClonalAdhesion, StaticAdhesion};
use crate::cellular_automata::CellularAutomata;
use crate::environment::{Environment, LatticeEntity};
use crate::io::io_manager::IoManager;
use crate::io::parameters::Parameters;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use std::error::Error;
use rand::distr::{Distribution, Uniform};
use crate::cell_container::CellContainer;
use crate::constants::NeighbourhoodType;
use crate::genome::{Genome, Grn};
use crate::io::movie_maker::MovieMaker;
use crate::positional::rect::Rect;
use crate::space::Space;

pub struct Model {
    pub env: Environment<NeighbourhoodType>,
    pub ca: CellularAutomata<ClonalAdhesion>,
    pub rng: Xoshiro256StarStar,
    pub io_manager: IoManager
}

// TODO! Implement new, which does not require Parameters
impl Model {
    pub fn run(&mut self, steps: u32) {
        log::info!("Starting simulation");
        for time_step in 0..=steps {
            let saved = self.io_manager.image_io(
                time_step,
                &self.env, 
                &self.ca.adhesion.clone_pairs
            );
            if let Err(e) = saved {
                log::warn!("Failed to save image at time step {time_step} with error `{e}`")
            }
            self.step(time_step);
        }
    }
    
    pub fn step(&mut self, time_step: u32) {
        self.ca.step(&mut self.env, &mut self.rng);
        if self.env.time_to_update(time_step) {
            self.env.cells.update_cells();
            let new_spins = self.env.reproduce();
            for spin in new_spins {
                self.ca.adhesion.update_clones(spin, &self.env);
                // We could also instead choose to mutate at a fix rate throughout the cell's life cycle
                if let LatticeEntity::SomeCell(cell) = self.env.cells.get_entity_mut(spin) {
                    cell.genome.attempt_mutate(&mut self.rng);
                } else { 
                    panic!("Newborn is not a cell")
                }
            }
        }
    }
}

impl TryFrom<Parameters> for Model {
    type Error = Box<dyn Error>;

    fn try_from(parameters: Parameters) -> Result<Self, Self::Error> {
        let rng = if parameters.general.seed == 0 {
            Xoshiro256StarStar::from_os_rng()
        } else {
            Xoshiro256StarStar::seed_from_u64(parameters.general.seed)
        };
        let mut model = Self {
            env: Environment::new(
                parameters.environment.update_period,
                parameters.environment.cell_search_radius,
                parameters.environment.max_cells,
                parameters.environment.enclose,
                CellContainer::new(
                    parameters.environment.cell.target_area,
                    parameters.environment.cell.div_area,
                    parameters.environment.cell.divide,
                    parameters.environment.cell.migrate,
                ),
                Space::new(
                    parameters.environment.width,
                    parameters.environment.height,
                ),
                NeighbourhoodType::new(parameters.environment.neigh_r)
            ),
            ca: CellularAutomata::new(
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
            ),
            io_manager: IoManager::new(
                &parameters.io.outdir,
                parameters.io.image_period,
                parameters.io.image_format.clone(),
                parameters.io.plots.clone(),
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
            ),
            rng
        };

        log::info!("Setting model up");
        log::info!("Creating output directory and copy of parameter file");
        model.io_manager.create_directories(parameters.io.replace_outdir)?;
        model.io_manager.create_parameters_file(&parameters)?;

        log::info!("Creating cells");
        let mut cell_count = 0;
        let cell_side = (parameters.environment.cell_start_area as f32).sqrt() as usize;
        for _ in 0..parameters.environment.starting_cells {
            let pos = model.env.space.cell_lattice.random_pos(&mut model.rng);
            let cell = model.env.spawn_rect_cell(
                Rect::new(
                    pos,
                    (pos.x + cell_side, pos.y + cell_side).into()
                ),
                // TODO!: Parameterize
                Grn::new(
                    1. / model.env.height() as f32,
                    2,
                    0.8,
                    0.8,
                    || Uniform::new(-1., 1.).unwrap().sample(&mut model.rng)
                )
            );
            if cell.is_some() {
                cell_count += 1;
            }
        }
        log::info!(
            "Created {} out of the {} cells requested", 
            cell_count, 
            parameters.environment.starting_cells
        );

        Ok(model)
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