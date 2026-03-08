## Cellulars: a fast and extensible implementation of the Cellular Potts Model

If you want full control over simulation code, this library is for you¹.
Due to Rust's incredibly flexible type system, modifying cell behavior using cellulars is easy and intuitive, 
while adding 0 overhead in execution time. 

Besides, cellulars was made to be really fast. 
It currently runs on a single thread and achieves high performance using algorithmic optimizations of the
Monte Carlo sampling method used traditionally for CPM.
In the near future, the main step function will be parallelized for an even larger speed-up (stay tuned).

Features currently include:
- Total control over the Hamiltonian functional (C).
- The Environment API makes it easy to design hybrid models that might require complex data structures (C).
- Plotting API to display information about the simulation (C).
- Real-time movie inspection of the simulation (C).
- Implementations of multiple boundary and neighborhood conditions (C).
- Write data about the simulation to use for analyses (C).
- Read simulation templates from image files (M).
- Robust backup system to restore simulations (M).

**C**: Features included in the cellulars lib.
**M**: Features included in the [model template](https://github.com/aleferna12/cellulars-template).

¹ There are also alternative implementations of the CPM for users who are not as comfortable in Rust and prefer an
out-of-the-box approach
(check out [Morpheus](https://morpheus.gitlab.io/) and [Artistoo](https://artistoo.net/), for example).

### Getting started:

The simplest CPM possible looks like this:

```rust
fn main() {
    // Initialize periodic boundary conditions
    let boundary = UnsafePeriodicBoundary::new(Rect::new(
        Pos::new(0., 0.),
        Pos::new(100., 100.)
    ));
    // Initialize an empty environment
    let mut env = Environment::new(
        cell_container![],
        Lattice::new(100, 100),
        MooreNeighborhood::new(1),
        Boundaries::new(boundary)
    );
    // Spawn a cell in a rectangular region of the environment
    let cell_rect = Rect::new(
        Pos::new(25, 25),
        Pos::new(75, 75)
    );
    env.spawn_cell(Cell::new_empty(cell_rect.area() as u32), cell_rect.iter_positions());

    // Initialize the Potts algorithm used to update the environment
    let potts = EdgePotts {
        boltz_t: 16.,
        size_lambda: 4.,
        adhesion: StaticAdhesion {
            cell_energy: 10.,
            medium_energy: 10.,
            solid_energy: 10.
        },
        bias: NoBias
    };
    let mut rng = ThreadRng::default();

    // Run 100k time steps of the simulation
    for _ in 0..100_000 {
        potts.step(&mut env, &mut rng);
    }
}
```

You will also need to attach a movie window to the simulation to see what is going on, see the full example at
[`cellulars/examples/basic.rs`](https://github.com/aleferna12/cellula-rs/blob/master/cellulars/examples/basic.rs).

For more examples, see the [examples](https://github.com/aleferna12/cellula-rs/blob/master/cellulars/examples) folder.
Also check out cellulars [documentation](https://docs.rs/cellulars/latest/cellulars/).

### Installation:

There are three ways to install cellulars, depending on the level of control you want over the simulation code.

#### 1: Clone the template repo (recommended for quick start)

To start with a template for a model that includes a CLI interface,
automatic IO managing, parameter parsing, and other quality of life features, head over to the
[model template page](https://github.com/aleferna12/cellulars-template) and follow installation instructions there.

#### 2: Add it as a dependency (recommended for experimentation/mixed project)

This will include only library code. Just add cellulars as a dependency to your cargo project with:

```commandline
cargo add cellulars
```

#### 3: Clone this repo (recommended for extreme granular control)

In principle, there should be few cases where this is useful, since the
lib was made highly extensible by design.
However, if for any reason you need to modify code from the cellulars library itself, you can clone this repo (which includes both the lib code
at `cellulars/` and the model code at `model/`):

```commandline
git clone https://github.com/aleferna12/cellula-rs <MY_PROJECT_NAME>
```
