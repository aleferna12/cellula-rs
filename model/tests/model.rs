use model::io::parameters::Parameters;
use model::io::parameters::PlotType;
use model::io::parameters::PlotType::*;
use model::model::Model;

fn make_test_model(plot: PlotType) -> anyhow::Result<Model> {
    let mut pars = Parameters::parse("examples/64_cells.toml")?;
    pars.io.outdir = format!("tests/{plot:?}").into();
    pars.io.plot.order = vec![Chem, Spin, plot, Border];
    pars.io.image_period = 50;
    pars.io.movie.show = false;
    pars.cell.update_period = 1;
    Model::initialise_from_parameters(pars)
}

#[test]
fn test_parse_model() -> anyhow::Result<()> {
    let model = make_test_model(Spin)?;
    assert!(model.pond.env.cells.n_cells() > 60);
    Ok(())
}

#[test]
fn test_run() {
    for plot in [Clones, CellType, Area, Center, ChemCenter] {
        let mut model = make_test_model(plot).unwrap();
        model.run_for(500);
    }
}