use std::cmp::min;
use std::hint::black_box;
use criterion::{criterion_group, criterion_main, Criterion};
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use evo_cpm::cell::Cell;
use evo_cpm::edge::Edge;
use evo_cpm::environment::Environment;
use evo_cpm::environment::LatticeEntity::*;
use evo_cpm::pos::Pos2D;

fn random_neighbour(env: &Environment, p: Pos2D<usize>, neigh_r: u8, rng: &mut impl Rng) -> Pos2D<usize> {
    let oldp = (p.x as i32, p.y as i32);
    let mut newp = oldp;
    let dist = neigh_r as i32;
    while oldp == newp {
        newp.0 = oldp.0 + rng.random_range(
            -min(dist, oldp.0)..min(dist + 1, env.cell_lattice.width() as i32 - oldp.0)
        );
        newp.1 = oldp.1 + rng.random_range(
            -min(dist, oldp.1)..min(dist + 1, env.cell_lattice.height() as i32 - oldp.1)
        );
    }
    Pos2D::new(newp.0 as usize, newp.1 as usize)
}

fn add_random_edge(env: &mut Environment, rng: &mut impl Rng) -> bool {
    let p1 = env.cell_lattice.random_pos(rng);
    let e = Edge::new(p1, random_neighbour(&env, p1, 1, rng));
    env.edge_book.insert(e)
}

fn replace_random_edges(n_edges: usize, env: &mut Environment, rng: &mut impl Rng) {
    for _ in 0..n_edges {
        let e1 = add_random_edge(env, rng);
        if e1 {
            let i = env.edge_book.random_index(rng);
            env.edge_book.remove_at(i);
        }
    }
}

fn bench_env(c: &mut Criterion) {
    let mut env = Environment::new(
        100, 
        100, 
        1, 
        0, 
        0, 
        0
    );
    let mut rng = Xoshiro256StarStar::seed_from_u64(1241254152);
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

    c.bench_function("lattice_entity_discriminant", |b| b.iter(|| {
        Medium::<&Cell>.spin();
        Solid::<&Cell>.spin();
    }));
}

criterion_group!(env_bench, bench_env);
criterion_main!(env_bench);
