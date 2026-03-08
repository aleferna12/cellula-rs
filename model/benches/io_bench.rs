use cellulars::traits::step::Step;
use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use model::io::parameters::{Parameters, PlotType};
use model::model::Model;
use strum::IntoEnumIterator;

fn make_model(params: Parameters) -> Model {
    let mut model = Model::new_from_parameters(params, None).unwrap();
    model.run_for(100);
    model
}

fn bench_io(c: &mut Criterion) {
    let mut params = Parameters::parse("config/1_cell.toml").unwrap();
    params.io.image_period = 1000000;
    #[cfg(feature = "movie-io")]
    if let Some(movie_params) = &mut params.io.movie {
        movie_params.show = false;
    }
    for plot in PlotType::iter() {
        c.bench_with_input(
            BenchmarkId::new("plot", format!("{plot:?}")),
            &plot,
            |b, i| {b.iter_batched_ref(
                || {
                    params.io.plot.order = vec![i.clone()].into();
                    make_model(params.clone())
                },
                |model| {
                    model.io.make_simulation_image(&model.pond.env);
                },
                BatchSize::LargeInput
            )}
        );
    }
}

criterion_group!(io_bench, bench_io);
criterion_main!(io_bench);