use std::hint::black_box;
use criterion::{criterion_group, Criterion};
use evo_cpm::lattice::Lattice;
use evo_cpm::pos::Edge;

fn add_random_edge(lat: &mut Lattice<u32>) -> bool {
    let p1 = lat.random_pos();
    let e = Edge::new(p1, lat.random_neighbour(&p1));
    if e.is_err() { return false }
    lat.insert_edge(e.unwrap())
}

fn replace_random_edges(n_edges: usize, lat: &mut Lattice<u32>) {
    for _ in 0..n_edges {
        let e = add_random_edge(lat);
        if e { lat.remove_random_edge(); }
    }
}

fn criterion_benchmark(c: &mut Criterion) {
    let mut lat = Lattice::new(20, 20);
    for _ in 0..lat.width * lat.height / 2 {
        add_random_edge(&mut lat);
    }
    c.bench_function("replace_edge", |b| { 
        b.iter(|| replace_random_edges(black_box(100_000), black_box(&mut lat)))
    });
}

criterion_group!(lattice, criterion_benchmark);