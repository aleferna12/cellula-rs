use criterion::BatchSize;
use criterion::{criterion_group, criterion_main, Criterion};
use evo_cpm::environment::LatticeEntity::*;
use evo_cpm::environment::{Environment, LatticeEntity};
use evo_cpm::positional::edge::Edge;
use evo_cpm::positional::pos::Pos;
use evo_cpm::positional::rect::Rect;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use std::cmp::min;
use std::hint::black_box;
use evo_cpm::positional::boundary::{Boundary, UnsafePeriodicBoundary};

fn random_neighbour(env: &Environment, p: Pos<usize>, neigh_r: u8, rng: &mut impl Rng) -> Pos<usize> {
    let oldp = (p.x as i32, p.y as i32);
    let mut newp = oldp;
    let dist = neigh_r as i32;
    while oldp == newp {
        newp.0 = oldp.0 + rng.random_range(
            -min(dist, oldp.0)..min(dist + 1, env.space.cell_lattice.width() as i32 - oldp.0)
        );
        newp.1 = oldp.1 + rng.random_range(
            -min(dist, oldp.1)..min(dist + 1, env.space.cell_lattice.height() as i32 - oldp.1)
        );
    }
    Pos::new(newp.0 as usize, newp.1 as usize)
}

fn add_random_edge(env: &mut Environment, rng: &mut impl Rng) -> bool {
    let p1 = env.space.cell_lattice.random_pos(rng);
    let e = Edge::new(p1, random_neighbour(env, p1, 1, rng));
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
    c.bench_function("lattice_entity_discriminant", |b| b.iter(|| {
        Medium.discriminant();
        Solid.discriminant();
    }));

    c.bench_function("replace_edges", |b| {
        b.iter_batched_ref(
            || {
                let mut env = Environment::new_empty_test(100, 100);
                let mut rng = Xoshiro256StarStar::seed_from_u64(1241254152);
                for _ in 0..env.space.cell_lattice.width() * env.space.cell_lattice.height() / 2 {
                    add_random_edge(&mut env, &mut rng);
                }
                (env, rng)
            },
            |(env, rng)| replace_random_edges(
                black_box(100_000),
                black_box(env),
                black_box(rng)
            ),
            BatchSize::SmallInput
        );
    });
    
    let pos_usize: [Pos<isize>; 2] = [Pos::new(20, 20), Pos::new(-20, -20)];
    let lat_bound = UnsafePeriodicBoundary::from(Rect::new((0, 0).into(), (40, 40).into()));
    c.bench_function("unsafe_periodic_boundary_usize", |b| {
        b.iter(
            || {
                lat_bound.valid_positions(pos_usize.into_iter())
            }
        );
    });

    let pos_usize: [Pos<f32>; 2] = [Pos::new(20., 20.), Pos::new(-20., -20.)];
    let lat_bound = UnsafePeriodicBoundary::from(Rect::new((0., 0.).into(), (40., 40.).into()));
    c.bench_function("unsafe_periodic_boundary_f32", |b| {
        b.iter(
            || {
                lat_bound.valid_positions(pos_usize.into_iter())
            }
        );
    });

    let mut env = Environment::new_empty_test(100, 100);
    env.spawn_rect_cell(Rect::new((10, 10).into(), (20, 20).into()));

    let mut group = c.benchmark_group("cell_positions");
    group.bench_function("contiguous_cell_positions", |b| {
        b.iter(|| {
            let cell = env.cells.get_entity(LatticeEntity::first_cell_spin()).unwrap_cell();
            assert_eq!(
                env.space.contiguous_cell_positions(cell, &env.neighbourhood).len(),
                cell.area as usize
            );
        })
    });
    group.bench_function("box_cell_positions", |b| {
        b.iter(|| {
            let cell = env.cells.get_entity(LatticeEntity::first_cell_spin()).unwrap_cell();
            assert_eq!(
                env.space.box_cell_positions(cell, 2.).len(),
                cell.area as usize
            );
        })
    });
}

criterion_group!(env_bench, bench_env);
criterion_main!(env_bench);
