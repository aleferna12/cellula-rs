mod env_bench;
mod model_bench;
mod io_bench;
use criterion::criterion_main;

criterion_main!(
    env_bench::env_bench,
    io_bench::io_bench,
    model_bench::model_1000mcs
);
