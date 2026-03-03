use cellulars::traits::step::Step;
use model::io::parameters::{Parameters, PlotType as PT};
use model::model::Model;

fn make_test_parameters() -> anyhow::Result<Parameters> {
    let mut params = Parameters::parse("examples/64_cells.toml")?;
    params.io.image_period = 64;
    params.io.data.cells_period = 512;
    params.io.data.lattice_period = 512;
    params.cell.update_period = 1;
    #[cfg(feature = "movie-io")]
    if let Some(movie_params) = &mut params.io.movie {
        movie_params.show = false;
    }
    Ok(params)
}

#[test]
fn test_plots() -> anyhow::Result<()> {
    for plot in [PT::CellType, PT::Area, PT::Center, PT::ChemCenter] {
        let mut params = make_test_parameters()?;
        params.io.outdir = format!("tests/out/plots/{plot:?}");
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
    params.io.outdir = "tests/out/templates/".to_string();

    let mut model = Model::new_from_parameters(params, Some("tests/fixtures/mig_div_templates.parquet".to_string()))?;
    model.run_for(512);
    Ok(())
}

#[test]
fn test_layout() -> anyhow::Result<()> {
    let mut params = make_test_parameters()?;
    params.io.outdir = "tests/out/layout/".to_string();

    let mut model = Model::new_from_layout(params, "tests/fixtures/squares_layout.png".to_string(), None)?;
    model.run_for(512);
    Ok(())
}

#[test]
fn test_layout_template() -> anyhow::Result<()> {
    let mut params = make_test_parameters()?;
    params.io.outdir = "tests/out/layout_template/".to_string();

    let mut model = Model::new_from_layout(
        params,
        "tests/fixtures/squares_layout.png".to_string(),
        Some("tests/fixtures/mig_div_templates.parquet".to_string())
    )?;
    model.run_for(512);
    Ok(())
}