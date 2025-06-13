mod env_bench;
mod model_bench;

use criterion::criterion_main;

criterion_main!(env_bench::env, model_bench::model);