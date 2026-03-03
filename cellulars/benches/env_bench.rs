use cellulars::base::cell::Cell;
use cellulars::base::environment::Environment;
use cellulars::constants::FloatType;
use cellulars::positional::boundaries::{Boundaries, Boundary, UnsafePeriodicBoundary};
use cellulars::positional::edge::Edge;
use cellulars::positional::neighborhood::MooreNeighborhood;
use cellulars::positional::pos::Pos;
use cellulars::positional::rect::Rect;
use cellulars::prelude::Cellular;
use cellulars::traits::habitable::Habitable;
use criterion::BatchSize;
use criterion::{criterion_group, criterion_main, Criterion};
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use std::cmp::min;
use std::default::Default;
use std::hint::black_box;

fn empty_env(width: FloatType, height: FloatType) -> Environment<Cell, MooreNeighborhood, UnsafePeriodicBoundary<FloatType>> {
    Environment::new_empty(
        MooreNeighborhood::new(1),
        Boundaries::new(UnsafePeriodicBoundary::new(Rect::new(
            (0., 0.).into(),
            (width, height).into()
        )))
    )
}

fn random_neighbor(
    env: &Environment<Cell, MooreNeighborhood, UnsafePeriodicBoundary<FloatType>>,
    p: Pos<usize>,
    neigh_r: u8,
    rng: &mut impl RngExt
) -> Pos<usize> {
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
    Pos::new(newp.0 as usize, newp.1 as usize)
}

fn add_random_edge(
    env: &mut Environment<Cell, MooreNeighborhood, UnsafePeriodicBoundary<FloatType>>,
    rng: &mut impl RngExt
) -> bool {
    let p1 = env.cell_lattice.random_pos(rng);
    let e = Edge::new(p1, random_neighbor(env, p1, 1, rng));
    env.edge_book.insert(e)
}

fn replace_random_edges(
    n_edges: usize,
    env: &mut Environment<Cell, MooreNeighborhood, UnsafePeriodicBoundary<FloatType>>,
    rng: &mut impl RngExt
) {
    for _ in 0..n_edges {
        let e1 = add_random_edge(env, rng);
        if e1 {
            let i = env.edge_book.random_index(rng);
            env.edge_book.remove_at(i);
        }
    }
}

fn bench_env(c: &mut Criterion) {
    c.bench_function("replace_edges", |b| {
        b.iter_batched_ref(
            || {
                let mut env = empty_env(100., 100.);
                let mut rng = Xoshiro256StarStar::seed_from_u64(1241254152);
                for _ in 0..env.cell_lattice.width() * env.cell_lattice.height() / 2 {
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
    let lat_bound = UnsafePeriodicBoundary::new(Rect::new((0, 0).into(), (40, 40).into()));
    c.bench_function("unsafe_periodic_boundary_usize", |b| {
        b.iter(
            || {
                lat_bound.valid_positions(pos_usize.into_iter())
            }
        );
    });

    let pos_usize: [Pos<FloatType>; 2] = [Pos::new(20., 20.), Pos::new(-20., -20.)];
    let lat_bound = UnsafePeriodicBoundary::new(Rect::new((0., 0.).into(), (40., 40.).into()));
    c.bench_function("unsafe_periodic_boundary_f32", |b| {
        b.iter(
            || {
                lat_bound.valid_positions(pos_usize.into_iter())
            }
        );
    });

    let mut env = empty_env(100., 100.);
    env.spawn_cell(
        Cell::new_empty(100),
        Rect::new((10, 10).into(), (20, 20).into()).iter_positions(),
    );

    let mut group = c.benchmark_group("cell_positions");
    group.bench_function("contiguous_cell_positions", |b| {
        b.iter(|| {
            let rel_cell = &env.cells[0];
            assert_eq!(
                env.search_cell_contiguous(rel_cell).len(),
                rel_cell.cell.area() as usize
            );
        })
    });
    group.bench_function("box_cell_positions", |b| {
        b.iter(|| {
            let rel_cell = &env.cells[0];
            assert_eq!(
                env.search_cell_box(rel_cell, 2.).len(),
                rel_cell.cell.area() as usize
            );
        })
    });
}

criterion_group!(env_bench, bench_env);
criterion_main!(env_bench);
