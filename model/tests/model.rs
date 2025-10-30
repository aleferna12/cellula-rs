use cellulars_lib::step::Step;
use model::io::parameters::Parameters;
use model::io::parameters::PlotType::*;
use model::model::Model;

#[test]
fn test_run() -> anyhow::Result<()> {
    for plot in [CellType, Area, Center, ChemCenter] {
        let mut params = Parameters::parse("examples/64_cells.toml")?;
        params.io.outdir = format!("tests/{plot:?}").into();
        params.io.plot.order = vec![Chem, Spin, plot, Border];
        params.io.image_period = 64;
        params.io.movie.show = false;
        params.cell.update_period = 1;
        let mut model = Model::initialise_from_parameters(params.clone())?;
        model.run_for(512);
        // For now we resort to lying abt the time to trick IoManager into writing info
        model.io.write_if_time(4096, &model.pond.env)?;
        
        let sim_dir = params.io.outdir.clone();
        params.io.outdir = params.io.outdir + "/resumed/";
        let mut res_model = Model::initialise_from_backup(params, sim_dir, 4096)?;
        res_model.run_for(128);
    }
    Ok(())
}