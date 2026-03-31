use cellulars_lib::step::Step;
use model::io::parameters::{Parameters, PlotType as PT};
use model::model::Model;


fn make_test_parameters() -> anyhow::Result<Parameters> {
    let mut params = Parameters::parse("examples/64_cells.toml")?;
    params.io.image_period = 64;
    params.io.data.cells_period = 512;
    params.io.data.cells_write_period = Some(params.io.data.cells_period);
    params.io.data.lattice_period = params.io.data.cells_period;
    params.io.data.act_period = params.io.data.cells_period;
    params.io.movie.show = false;
    Ok(params)
}

#[test]
fn test_run() -> anyhow::Result<()> {
    for plot in [PT::Act, PT::Center, PT::ChemCenter, PT::Area, PT::RelChem] {
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

    let mut model = Model::new_from_parameters(params, Some("tests/big_small_templates.parquet".to_string()))?;
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
        Some("tests/big_small_templates.parquet".to_string())
    )?;
    model.run_for(512);
    Ok(())
}