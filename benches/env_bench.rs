use std::hint::black_box;
use criterion::{criterion_group, Criterion};
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use evo_cpm::environment::Environment;
use evo_cpm::pos::Edge;
use evo_cpm::utils::TEST_SEED;

fn add_random_edge(env: &mut Environment, rng: &mut impl Rng) -> bool {
    let p1 = env.cell_lattice.random_pos(rng);
    let e = Edge::new(p1, env.random_neighbour(&p1, 1, rng), 1);
    if e.is_err() { return false }
    env.insert_edge(e.unwrap())
}

fn replace_random_edges(n_edges: usize, env: &mut Environment, rng: &mut impl Rng) {
    for _ in 0..n_edges {
        let e = add_random_edge(env, rng);
        if e { env.remove_random_edge(rng); }
    }
}

fn random_edges(c: &mut Criterion) {
    let mut env = Environment::new(100, 100, 1);
    let mut rng = Xoshiro256StarStar::seed_from_u64(TEST_SEED);
    for _ in 0..env.cell_lattice.width() * env.cell_lattice.height() / 2 {
        add_random_edge(&mut env, &mut rng);
    }
    c.bench_function("replace_edges", |b| {
        b.iter(|| replace_random_edges(
            black_box(100_000),
            black_box(&mut env),
            black_box(&mut rng))
        )
    });
}

criterion_group!(env, random_edges);