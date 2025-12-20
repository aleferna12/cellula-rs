use cellulars_lib::step::Step;
use model::io::parameters::{Parameters, PlotType as PT};
use model::model::Model;

#[test]
fn test_run() -> anyhow::Result<()> {
    for plot in [PT::CellType, PT::Area, PT::Center, PT::ChemCenter] {
        let mut params = Parameters::parse("examples/64_cells.toml")?;
        params.io.outdir = format!("tests/{plot:?}");
        params.io.plot.order = vec![PT::Chem, PT::Spin, plot, PT::Border].into();
        params.io.image_period = 64;
        #[cfg(feature = "movie")]
        if let Some(movie_params) = &mut params.io.movie {
            movie_params.show = false;
        }
        params.cell.update_period = 1;
        let mut model = Model::new_from_parameters(params.clone())?;
        model.run_for(512);
        // For now we resort to lying abt the time to trick IoManager into writing info
        model.io.write_if_time(4096, &model.pond.env)?;
        
        let sim_dir = params.io.outdir.clone();
        params.io.outdir += "/resumed/";
        let mut res_model = Model::new_from_backup(params, sim_dir, 4096)?;
        res_model.run_for(128);
    }
    Ok(())
}