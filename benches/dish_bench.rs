use std::hint::black_box;
use criterion::{criterion_group, Criterion};
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use evo_cpm::dish::Dish;
use evo_cpm::pos::Edge;

fn add_random_edge(dish: &mut Dish, rng: &mut impl Rng) -> bool {
    let p1 = dish.cell_lattice.random_pos();
    let e = Edge::new(p1, dish.random_neighbour(&p1, 1, rng), 1);
    if e.is_err() { return false }
    dish.insert_edge(e.unwrap())
}

fn replace_random_edges(n_edges: usize, dish: &mut Dish, rng: &mut impl Rng) {
    for _ in 0..n_edges {
        let e = add_random_edge(dish, rng);
        if e { dish.remove_random_edge(rng); }
    }
}

fn random_edges(c: &mut Criterion) {
    let mut dish = Dish::new(100, 100, 1);
    let mut rng = Xoshiro256StarStar::from_os_rng();
    for _ in 0..dish.cell_lattice.width * dish.cell_lattice.height / 2 {
        add_random_edge(&mut dish, &mut rng);
    }
    c.bench_function("replace_edges", |b| {
        b.iter(|| replace_random_edges(
            black_box(100_000),
            black_box(&mut dish),
            black_box(&mut rng))
        )
    });
}

criterion_group!(dish, random_edges);