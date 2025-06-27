use std::cmp::max;
use criterion::{Criterion, criterion_group, criterion_main, BenchmarkId, BatchSize};
use evo_cpm::io::read_config;
use evo_cpm::model::Model;
use std::hint::black_box;
use std::time::Duration;
use evo_cpm::parameters::Parameters;

/// Builds all example models.
fn get_example_models() -> Vec<(String, Parameters)> {
    std::fs::read_dir("examples")
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
                    read_config(p).unwrap(),
                ))
            }
            _ => None,
        })
        .collect()
}

fn bench_examples(c: &mut Criterion, time_steps: u32) {
    for (example, parameters) in get_example_models() {
        c.bench_with_input(
            BenchmarkId::new("examples", format!("{}/{}mcs", example, time_steps)),
            &parameters,
            |b, parameters| {
                b.iter_batched_ref(
                    || {
                        let mut params = parameters.clone();
                        params.outdir = format!("benches/model_output/{}", params.outdir);
                        // Ensures that a single image will be saved, 
                        // either after the setup run or the whole simulation
                        params.image_period = max(time_steps, 100);
                        let mut model = Model::new(params);
                        model.setup().unwrap();
                        model.run(100);
                        model
                    },
                    |model| {
                        model.run(black_box(time_steps))
                    },
                    BatchSize::SmallInput
                )
            }
        );
    }
}

fn bench_examples_1000mcs(c: &mut Criterion) {
    bench_examples(c, 1000);
}

fn bench_examples_1mcs(c: &mut Criterion) {
    bench_examples(c, 1);
}

criterion_group!(
    name = model_1000mcs;
    config = Criterion::default()
        .sample_size(10)
        .measurement_time(Duration::from_secs(10));
    targets = bench_examples_1000mcs
);

criterion_group!(model_1mcs, bench_examples_1mcs);

criterion_main!(model_1mcs, model_1000mcs);
