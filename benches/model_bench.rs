use std::cmp::max;
use criterion::{Criterion, criterion_group, criterion_main, BenchmarkId, BatchSize};
use evo_cpm::model::Model;
use std::hint::black_box;
use std::path::Path;
use std::time::Duration;
use evo_cpm::parameters::Parameters;

/// Builds all example models.
fn find_parameters(parent_dir: impl AsRef<Path>) -> Vec<(String, Parameters)> {
    std::fs::read_dir(parent_dir)
        .unwrap()
        .filter_map(|entry| match entry {
            Ok(e) => {
                let p = e.path();
                if !p.is_file() || p.extension().unwrap().to_ascii_lowercase() != "toml" {
                    return None;
                }
                let bench_name = p.file_stem()
                    .unwrap()
                    .to_str()
                    .unwrap();
                Some((
                    bench_name.to_string(),
                    Parameters::from_file(p).unwrap(),
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
            BenchmarkId::new(function_prefix, format!("{}/{}mcs", file_name, time_steps)),
            &parameters,
            |b, parameters| {
                b.iter_batched_ref(
                    || {
                        let mut params = parameters.clone();
                        params.io.outdir = format!("benches/model_outputs/{}", params.io.outdir);
                        // Ensures that a single image will be saved, 
                        // either after the setup run or the whole simulation
                        params.io.image_period = max(time_steps, 100);
                        let mut model = Model::try_from(params).unwrap();
                        model.run(100);
                        model
                    },
                    |model| {
                        model.run(black_box(time_steps))
                    },
                    BatchSize::SmallInput,
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

criterion_group!(
    name = model_1000mcs;
    config = Criterion::default()
        .sample_size(10)
        .measurement_time(Duration::from_secs(10));
    targets = bench_param_files_1000mcs
);

criterion_group!(model_1mcs, bench_param_files_1mcs);

criterion_main!(model_1mcs, model_1000mcs);
