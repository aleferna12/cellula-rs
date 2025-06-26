// TODO!: fix (should read from the example files)
use std::error::Error;
use std::hint::black_box;
use std::time::Duration;
use clap::Parser;
use criterion::{criterion_group, criterion_main, Criterion};
use evo_cpm::model::Model;
use evo_cpm::parameters::Parameters;

/// This model should override all parameters that can have an effect on performance 
/// (do not depend on changeable defaults).
///
/// In the future, we can implement a test_config.toml with all of these.
pub fn test_parameters() -> Parameters {
    let mut params = Parameters::parse_from([""]);
    params.seed = 123451;
    params.replace_outdir = true;
    params.width = 1000;
    params.height = 1000;
    params.cell_start_area = 50;
    params.cell_target_area = 50;
    params.neigh_r = 1;
    params.boltz_t = 16.;
    params.size_lambda = 1.;
    params.cell_energy = 16.;
    params.medium_energy = 16.;
    params.solid_energy = 16.;
    params
}

fn hundred_cells_model() -> Result<Model, Box<dyn Error>> {
    let mut model = Model::new(test_parameters());
    model.setup()?;
    model.run(1000)?;
    Ok(model)
} 

fn single_cell_model() -> Result<Model, Box<dyn Error>> {
    let mut params = test_parameters();
    params.n_cells = 1;
    let mut model = Model::new(params);
    model.setup()?;
    Ok(model)
}

fn bench_model_1000mcs(c: &mut Criterion) {
    let mut sc_model = single_cell_model().unwrap();
    c.bench_function("single_cell_1000mcs", |b| {
        b.iter(|| {
            sc_model.run(black_box(1_000)).unwrap()
        })
    });

    let mut f_model = hundred_cells_model().unwrap();
    c.bench_function("full_model_1000mcs", |b| {
        b.iter(|| {
            f_model.run(black_box(1_000)).unwrap()
        })
    });
}

fn bench_model_1mcs(c: &mut Criterion) {
    let mut sc_model = single_cell_model().unwrap();
    c.bench_function("single_cell_1mcs", |b| {
        b.iter(|| {
            sc_model.run(black_box(1)).unwrap()
        })
    });

    let mut f_model = hundred_cells_model().unwrap();
    c.bench_function("full_model_1mcs", |b| {
        b.iter(|| {
            f_model.run(black_box(1)).unwrap()
        })
    });
}

criterion_group!(
    name = model_1000mcs;
    config = Criterion::default()
        .sample_size(10)
        .measurement_time(Duration::from_secs(10));
    targets = bench_model_1000mcs
);

criterion_group!(model_1mcs, bench_model_1mcs);

criterion_main!(model_1mcs, model_1000mcs);
