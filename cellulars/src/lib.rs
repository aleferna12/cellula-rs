/*!
A very fast, extensible CPM implementation.

## Getting started:

A minimal CPM implementation can be found in the
[basics example](https://github.com/aleferna12/cellula-rs/blob/master/cellulars/examples/basic.rs),
which should give you an idea of how to get started.

Also check out the other [examples](https://github.com/aleferna12/cellula-rs/blob/master/cellulars/examples).

Finally, we provide a [template repo](https://github.com/aleferna12/cellulars-template)
with additional features that do not fit in the library to reduce boilerplate.
These include a command line tool, a parameter parser and automatic data-write/backup.
*/

pub mod lattice;
pub mod positional;
pub mod static_adhesion;
pub mod cell_container;
pub mod symmetric_table;
pub mod constants;
pub mod spin;
pub mod traits;
pub mod base;
#[cfg(any(feature = "data-io", feature = "image-io"))]
pub mod io;
pub mod prelude;
pub mod empty_cell;