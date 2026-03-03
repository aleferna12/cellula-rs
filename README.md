## Cellulars: a fast and extensible implementation of the Cellular Potts Model

If you want full control over simulation code, this library is for you*.
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

\* there are also alternative implementations of the CPM for users that are not as comfortable in Rust and prefer an
out-of-the-box approach
(check out [Morpheus](https://morpheus.gitlab.io/) and [Artistoo](https://artistoo.net/), for example).

### Installation:

There are three ways to install cellulars, depending on the level of control you want over the simulation code.

#### 1: Clone the template repo (recommended for full CPM project)

To start with a template for a model that includes a CLI interface,
automatic IO managing and other quality of life features, clone the repo at 
(https://github.com/aleferna12/cellulars-template) and start from there (read the README).

#### 2: Add it as a dependency (recommended for experimentation/mixed project)

This will include only library code (so no fancy built-in IO/CLI features, which you will need to implement yourself).
Just add cellulars as a dependency to your cargo project with:

```commandline
cargo add cellulars
```

#### 3: Clone this repo (recommended for extreme granular control)

If you need to modify code from the cellulars library itself, you can clone this repo (which includes both the lib code
at `cellulars/` and the model code at `model/`). In principle there should be few cases where this is useful, since the
lib was made highly extensible by design.