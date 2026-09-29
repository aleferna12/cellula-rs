use cellulars_lib::traits::cellular::{Alive, Cellular};
use cellulars_lib::traits::step::Step;
use cellulars_lib::traits::track_perimeter::{measure_perimeters, TrackPerimeter};
use model::io::parameters::{Parameters, PlotType as PT};
use model::model::Model;

fn make_test_parameters() -> anyhow::Result<Parameters> {
    let mut params = Parameters::parse("examples/64_cells.toml")?;
    params.io.image_period = 64;
    params.io.data.cells_period = 512;
    params.io.data.lattice_period = 512;
    params.io.kinect = None;
    params.cell.update_period = 1;
    #[cfg(feature = "movie")]
    if let Some(movie_params) = &mut params.io.movie {
        movie_params.show = false;
    }
    Ok(params)
}

#[test]
fn test_run() -> anyhow::Result<()> {
    for plot in [PT::CellType, PT::Area, PT::Center, PT::ChemCenter] {
        let mut params = make_test_parameters()?;
        params.io.outdir = format!("tests/plots/{plot:?}");
        params.io.plot.order = vec![PT::Chem, PT::Spin, plot, PT::Border].into();

        let mut model = Model::new_from_parameters(params.clone(), None)?;
        model.run_for(513);
        
        let sim_dir = params.io.outdir.clone();
        params.io.outdir += "/resumed/";
        let mut res_model = Model::new_from_backup(params, sim_dir, 512)?;
        res_model.run_for(128);
    }
    Ok(())
}

#[test]
fn test_templates() -> anyhow::Result<()> {
    let mut params = make_test_parameters()?;
    params.io.outdir = "tests/templates/".to_string();

    let mut model = Model::new_from_parameters(params, Some("tests/mig_div_templates.parquet".to_string()))?;
    model.run_for(512);
    Ok(())
}

#[test]
fn test_layout() -> anyhow::Result<()> {
    let mut params = make_test_parameters()?;
    params.io.outdir = "tests/layout/".to_string();

    let mut model = Model::new_from_layout(params, "tests/squares_layout.png".to_string(), None)?;
    model.run_for(512);
    Ok(())
}

#[test]
fn test_layout_template() -> anyhow::Result<()> {
    let mut params = make_test_parameters()?;
    params.io.outdir = "tests/layout_template/".to_string();

    let mut model = Model::new_from_layout(
        params,
        "tests/squares_layout.png".to_string(),
        Some("tests/mig_div_templates.parquet".to_string())
    )?;
    model.run_for(512);
    Ok(())
}
/// The perimeters tracked by the cells must match the ones measured from the lattice,
/// otherwise the perimeter constraint is penalising the wrong deviations.
#[test]
fn test_perimeter_tracking() -> anyhow::Result<()> {
    let mut params = make_test_parameters()?;
    params.io.outdir = "tests/out/perimeter".into();
    params.cell.target_perimeter = 170;
    params.potts.perimeter_lambda = 0.5;

    let mut model = Model::new_from_parameters(params, None)?;
    model.run_for(256);

    let env = model.my_pond.env();
    let measured = measure_perimeters(&env.env);
    for rel_cell in env.env.cells.iter() {
        assert_eq!(
            rel_cell.cell.perimeter(),
            measured[rel_cell.index as usize],
            "perimeter of cell {} drifted from the lattice",
            rel_cell.index
        );
    }
    Ok(())
}

/// Cells must come out of initialisation with a unit migration direction and a staggered persistence
/// clock, and must keep their directions normalised as they turn.
#[test]
fn test_persistent_migration() -> anyhow::Result<()> {
    let mut params = make_test_parameters()?;
    params.io.outdir = "tests/out/persistence".into();
    params.cell.persistence_duration = 8;
    params.potts.persistence_mu = 2.;

    let mut model = Model::new_from_parameters(params, None)?;
    for rel_cell in model.my_pond.env().env.cells.iter() {
        assert!(rel_cell.cell.has_direction(), "cell {} was not given a direction", rel_cell.index);
        assert!(rel_cell.cell.persistence_time() < rel_cell.cell.persistence_duration);
    }
    // Cells start with a staggered clock, so they should not all be turning on the same time-step
    let clocks: Vec<_> = model
        .my_pond
        .env()
        .env
        .cells
        .iter()
        .map(|rel_cell| rel_cell.cell.persistence_time())
        .collect();
    assert!(clocks.iter().any(|time| time != &clocks[0]), "persistence clocks were not staggered");

    model.run_for(256);
    for rel_cell in model.my_pond.env().env.cells.iter() {
        if !rel_cell.cell.is_alive() {
            continue;
        }
        let target_vec = rel_cell.cell.target_vec();
        let hyp = target_vec.x.hypot(target_vec.y);
        assert!(
            (hyp - 1.).abs() < 1e-5,
            "direction of cell {} is not a unit vector (norm was {hyp})",
            rel_cell.index
        );
    }
    Ok(())
}

/// Migration must displace cells further than diffusion alone would.
#[test]
fn test_persistence_displaces_cells() -> anyhow::Result<()> {
    let run = |persistence_mu: f64, outdir: &str| -> anyhow::Result<f64> {
        let mut params = make_test_parameters()?;
        params.io.outdir = outdir.into();
        params.cell.divide = false;
        // Isolate persistent migration from the chemical gradient
        params.potts.chemotaxis_mu = 0.;
        params.cell.persistence_duration = 32;
        params.potts.persistence_mu = persistence_mu as _;

        let mut model = Model::new_from_parameters(params, None)?;
        let starts: Vec<_> = model
            .my_pond
            .env()
            .env
            .cells
            .iter()
            .map(|rel_cell| rel_cell.cell.center())
            .collect();
        model.run_for(512);

        let env = model.my_pond.env();
        let mut total = 0.;
        let mut counted = 0;
        for (rel_cell, start) in env.env.cells.iter().zip(starts) {
            if !rel_cell.cell.is_alive() {
                continue;
            }
            let (dx, dy) = cellulars_lib::positional::boundaries::Boundary::displacement(
                &env.env.bounds.boundary,
                start,
                rel_cell.cell.center()
            );
            total += dx.hypot(dy) as f64;
            counted += 1;
        }
        Ok(total / counted as f64)
    };

    let still = run(0., "tests/out/persistence_off")?;
    let migrating = run(4., "tests/out/persistence_on")?;
    assert!(
        migrating > still,
        "persistent migration did not displace cells further than no migration \
         (migrating travelled {migrating:.2}, still travelled {still:.2})"
    );
    Ok(())
}
