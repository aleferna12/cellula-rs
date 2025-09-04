mod model_bench;
mod io_bench;
use criterion::criterion_main;

criterion_main!(
    io_bench::io_bench,
    model_bench::model_1000mcs
);
