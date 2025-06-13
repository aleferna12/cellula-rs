use std::hint::black_box;
use std::time::Duration;
use clap::Parser;
use criterion::{criterion_group, Criterion};
use evo_cpm::model::Model;
use evo_cpm::parameters::Parameters;
use evo_cpm::utils::TEST_SEED;

fn run_single_cell(steps: u32) {
    let mut model = Model::new(Parameters::parse_from(["", "--seed", &TEST_SEED.to_string()]));
    model.setup();
    model.run(steps);
}

fn bench_single_cell(c: &mut Criterion) {
    c.bench_function("single_cell_1e4", |b| {
        b.iter(|| {
            run_single_cell(black_box(10_000))
        })
    });
}

fn config() -> Criterion {
    Criterion::default()
        .sample_size(10)
        .measurement_time(Duration::from_secs(10))
}

criterion_group!(
    name = model;
    config = config();
    targets = bench_single_cell
);