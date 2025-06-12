mod env_bench;

use criterion::criterion_main;

criterion_main!(env_bench::env);