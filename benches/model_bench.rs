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
fn full_model() -> Model {
    let mut params = Parameters::parse_from([""]);
    params.seed = 123451;
    params.n_cells = 100;
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
    Model::new(params)
}

fn run_full_model(steps: u32) {
    let mut model = full_model();
    model.setup();
    model.run(1000);
    model.run(steps);
} 

fn run_single_cell(steps: u32) {
    let mut model = full_model();
    model.parameters.n_cells = 1;
    model.setup();
    model.run(1000);
    model.run(steps);
}

fn bench_model(c: &mut Criterion) {
    c.bench_function("single_cell_1000", |b| {
        b.iter(|| {
            run_single_cell(black_box(1_000))
        })
    });
    c.bench_function("full_model_1000", |b| {
        b.iter(|| {
            run_full_model(black_box(1_000))
        })
    });
}

fn config() -> Criterion {
    Criterion::default()
        .sample_size(10)
        .measurement_time(Duration::from_secs(10))
}

criterion_group!(
    name = model_bench;
    config = config();
    targets = bench_model
);
criterion_main!(model_bench);
