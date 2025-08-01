use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use evo_cpm::io::parameters::Parameters;
use evo_cpm::model::Model;
use std::cmp::max;
use std::hint::black_box;
use std::path::Path;
use std::time::Duration;

/// Builds all example models.
fn find_parameters(parent_dir: impl AsRef<Path>) -> Vec<(String, Parameters)> {
    std::fs::read_dir(parent_dir)
        .unwrap()
        .filter_map(|entry| match entry {
            Ok(e) => {
                let p = e.path();
                if !p.is_file() || !p.extension().unwrap().eq_ignore_ascii_case("toml") {
                    return None;
                }
                let bench_name = p.file_stem()
                    .unwrap()
                    .to_str()
                    .unwrap();
                Some((
                    bench_name.to_string(),
                    Parameters::parse(p).unwrap(),
                ))
            }
            _ => None,
        })
        .collect()
}

/// Searches `parent_dir` for config files, reads them, 
/// and them for each of them makes a benchmark called `function_prefix`/`file_name`/`time_steps`mcs.
fn bench_param_files(
    c: &mut Criterion,
    function_prefix: &str,
    parent_dir: impl AsRef<Path>,
    time_steps: u32
) {
    let parent_dir = parent_dir.as_ref();
    for (file_name, parameters) in find_parameters(parent_dir) {
        c.bench_with_input(
            BenchmarkId::new(function_prefix, format!("{file_name}/{time_steps}mcs")),
            &parameters,
            |b, parameters| {
                b.iter_batched_ref(
                    || {
                        let mut params = parameters.clone();
                        // Ensures that a single image will be saved, 
                        // either after the setup run or the whole simulation
                        params.io.image_period = max(time_steps, 100);
                        params.io.movie.show = false;
                        let mut model = Model::initialise_from_parameters(params).unwrap();
                        model.run_for(100);
                        model
                    },
                    |model| {
                        model.run_for(black_box(time_steps))
                    },
                    BatchSize::LargeInput,
                )
            },
        );
    }
}

fn bench_param_files_1000mcs(c: &mut Criterion) {
    bench_param_files(c, "examples", "./examples", 1000);
    bench_param_files(c, "models", "./benches/model_files", 1000);
}

fn bench_param_files_1mcs(c: &mut Criterion) {
    bench_param_files(c, "examples", "./examples", 1);
    bench_param_files(c, "models", "./benches/model_files", 1);
}

fn bench_slow(c: &mut Criterion) {
    c.bench_function("clonal_dev", |b| {
        b.iter_batched_ref(
            || {
                let mut params = Parameters::parse("examples/1_cell.toml").unwrap();
                params.io.image_period = 50_000;
                params.io.movie.show = false;
                let mut model = Model::initialise_from_parameters(params).unwrap();
                model.run_for(50_000);
                model
            },
            |model| {
                model.run_for(black_box(1_000));
            },
            BatchSize::LargeInput
        )
    });

    c.bench_function("large_lattice/10000mcs", |b| {
        b.iter_batched_ref(
            || {
                let mut params = Parameters::parse(
                    "./benches/model_files/large_lattice.toml"
                ).unwrap();
                params.io.image_period = 10_000;
                params.io.movie.show = false;
                let mut model = Model::initialise_from_parameters(params).unwrap();
                model.run_for(100);
                model
            },
            |model| {
                model.run_for(black_box(10_000));
            },
            BatchSize::LargeInput
        )
    });
}

criterion_group!(
    name = model_1000mcs;
    config = Criterion::default()
        .sample_size(10)
        .measurement_time(Duration::from_secs(15));
    targets = bench_param_files_1000mcs
);

criterion_group!(model_1mcs, bench_param_files_1mcs);

criterion_group!(
    name = slow;
    config = Criterion::default()
        .sample_size(10)
        .measurement_time(Duration::from_secs(240));
    targets = bench_slow
);

criterion_main!(model_1mcs, model_1000mcs, slow);
