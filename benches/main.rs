mod env_bench;
mod model_bench;
use criterion::criterion_main;

criterion_main!(env_bench::env_bench, model_bench::model_bench);
