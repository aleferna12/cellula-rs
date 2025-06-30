use std::hint::black_box;
use criterion::{criterion_group, criterion_main, Criterion};
use evo_cpm::ca::{AdhesionManager, StaticAdhesion};
use evo_cpm::cell::Cell;
use evo_cpm::environment::LatticeEntity;

fn bench_ca(c: &mut Criterion) {
    // These don't seem to be very reliable
    let mut group = c.benchmark_group("adhesion");
    let cell = Cell::new(10, 10);
    let some_cell = LatticeEntity::SomeCell(&cell);
    
    let unboxed_adh = StaticAdhesion { cell_energy: 16., medium_energy: 16., solid_energy: 16. };
    group.bench_function("unboxed_adhesion", |b| {
        b.iter(|| {
            unboxed_adh.adhesion_energy(black_box(some_cell), black_box(some_cell));
            unboxed_adh.adhesion_energy(black_box(some_cell), black_box(LatticeEntity::Medium));
            unboxed_adh.adhesion_energy(black_box(some_cell), black_box(LatticeEntity::Solid));
            unboxed_adh.adhesion_energy(black_box(LatticeEntity::Medium), black_box(LatticeEntity::Medium));
            unboxed_adh.adhesion_energy(black_box(LatticeEntity::Medium), black_box(LatticeEntity::Solid));
        })
    });

    let boxed_adh = Box::new(unboxed_adh);
    group.bench_function("boxed_adhesion", |b| {
        b.iter(|| {
            boxed_adh.adhesion_energy(black_box(some_cell), black_box(some_cell));
            boxed_adh.adhesion_energy(black_box(some_cell), black_box(LatticeEntity::Medium));
            boxed_adh.adhesion_energy(black_box(some_cell), black_box(LatticeEntity::Solid));
            boxed_adh.adhesion_energy(black_box(LatticeEntity::Medium), black_box(LatticeEntity::Medium));
            boxed_adh.adhesion_energy(black_box(LatticeEntity::Medium), black_box(LatticeEntity::Solid));
        })
    });
}

criterion_group!(ca_bench, bench_ca);
criterion_main!(ca_bench);
