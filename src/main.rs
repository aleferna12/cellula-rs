use crate::pos::{Edge};

mod lattice;
mod pos;

fn main() {
    let l = lattice::Lattice::<u32>::new(20, 20);
    println!("{}", l[pos::Pos2D::new(3, 3)]);
    println!("{}", l[(3, 3)]);
}
