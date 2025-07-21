use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use evo_cpm::io::parameters::{Parameters, PlotType};
use evo_cpm::model::Model;
use strum::IntoEnumIterator;

fn make_model() -> Model {
    let mut params = Parameters::parse("examples/64_cells.toml").unwrap();
    params.io.image_period = 1000000;
    params.io.movie.show = false;
    let mut model = Model::try_from(params).unwrap();
    model.run(100);
    model
}

fn bench_io(c: &mut Criterion) {
    for plot in PlotType::iter() {
        c.bench_with_input(
            BenchmarkId::new("plot", format!("{plot:?}")),
            &plot,
            |b, i| {b.iter_batched_ref(
                || {
                    let mut model = make_model();
                    model.io_manager.plots.order = vec![i.clone()];
                    model
                },
                |model| {
                    model.io_manager.simulation_image(&model.env, &model.ca.adhesion.clone_pairs).unwrap();
                },
                BatchSize::LargeInput
            )}
        );
    }
}

criterion_group!(io_bench, bench_io);
criterion_main!(io_bench);