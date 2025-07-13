mod env_bench;
mod model_bench;
mod ca_bench;
use criterion::criterion_main;

criterion_main!(
    ca_bench::ca_bench,
    env_bench::env_bench,
    model_bench::model_1000mcs
);
