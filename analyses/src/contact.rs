use anyhow::Context;
use cellulars_lib::positional::boundaries::{Boundary, FixedBoundary};
use cellulars_lib::positional::neighbourhood::{MooreNeighbourhood, Neighbourhood};
use cellulars_lib::positional::rect::Rect;
use cellulars_lib::spin::Spin;
use model::io::io_manager::IoManager;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use cellulars_lib::lattice::Lattice;
use cellulars_lib::positional::pos::Pos;

#[pymodule]
pub mod contact {
    #[pymodule_export]
    use super::local_act;
    #[pymodule_export]
    use super::geom_act;
}

#[pyfunction]
pub fn local_act(
    cell_lattice_path: &str,
    act_lattice_path: &str,
    width: usize,
    height: usize
) -> PyResult<ActVecs> {
    let (clat, alat) = clat_alat(
        cell_lattice_path,
        act_lattice_path,
        width,
        height
    ).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let neighborhood = MooreNeighbourhood::new(1);
    let bound = FixedBoundary::new(Rect::new(
        clat.rect.min.to_isize(),
        clat.rect.max.to_isize()
    ));
    let mut act_vecs = ActVecs::default();
    for pos in clat.iter_positions() {
        let spin = clat[pos];
        let Spin::Some(cell_index) = spin else {
            continue;
        };
        for neigh in neighborhood.neighbours(pos.to_isize()) {
            let Some(valid_neigh) = bound.valid_pos(neigh) else {
                continue;
            };
            let neigh_spin = clat[valid_neigh.to_usize()];
            if let Spin::Some(neighbor_index) = neigh_spin && cell_index == neighbor_index {
                continue;
            }
            let vec = act_vecs.vec_ref_mut(neigh_spin);
            vec.push(alat[pos] as f64);
        }
    }

    Ok(act_vecs)
}

#[pyfunction]
pub fn geom_act(
    cell_lattice_path: &str,
    act_lattice_path: &str,
    width: usize,
    height: usize
) -> PyResult<ActVecs> {
    let (clat, alat) = clat_alat(
        cell_lattice_path,
        act_lattice_path,
        width,
        height
    ).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let neighborhood = MooreNeighbourhood::new(1);
    let bound = FixedBoundary::new(Rect::new(
        clat.rect.min.to_isize(),
        clat.rect.max.to_isize()
    ));
    let mut act_vecs = ActVecs::default();
    for pos in clat.iter_positions() {
        let spin = clat[pos];
        let Spin::Some(cell_index) = spin else {
            continue;
        };
        let act = geom_mean_act(pos, &clat, &alat, &neighborhood, &bound);
        for neigh in neighborhood.neighbours(pos.to_isize()) {
            let Some(valid_neigh) = bound.valid_pos(neigh) else {
                continue;
            };
            let lat_neigh = valid_neigh.to_usize();
            let neigh_spin = clat[lat_neigh];
            if let Spin::Some(neighbor_index) = neigh_spin && cell_index == neighbor_index {
                continue;
            }
            let vec = act_vecs.vec_ref_mut(neigh_spin);
            let neigh_act = geom_mean_act(lat_neigh, &clat, &alat, &neighborhood, &bound);
            vec.push(neigh_act - act);
        }
    }

    Ok(act_vecs)
}

fn geom_mean_act(
    pos: Pos<usize>,
    clat: &Lattice<Spin>,
    alat: &Lattice<u32>,
    neighborhood: &impl Neighbourhood,
    bound: &impl Boundary<Coord = isize>
) -> f64 {
    let spin = clat[pos];
    let mut neigh_count = 0;
    let mut prod = 1;
    for neigh in neighborhood.neighbours(pos.to_isize()) {
        let Some(valid_neigh) = bound.valid_pos(neigh) else {
            continue;
        };
        let lat_neigh = valid_neigh.to_usize();
        if spin != clat[lat_neigh] {
            continue;
        }
        prod *= alat[lat_neigh];
        neigh_count += 1;
    }
    (prod as f64).powf(1. / neigh_count as f64)
}

#[pyclass]
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ActVecs {
    #[pyo3(get, set)]
    pub cell: Vec<f64>,
    #[pyo3(get, set)]
    pub medium: Vec<f64>,
    #[pyo3(get, set)]
    pub solid: Vec<f64>
}

impl ActVecs {
    #[allow(dead_code)]
    fn vec_ref(&self, spin: Spin) -> &Vec<f64> {
        match spin {
            Spin::Medium => &self.medium,
            Spin::Solid => &self.solid,
            Spin::Some(_) => &self.cell
        }
    }

    fn vec_ref_mut(&mut self, spin: Spin) -> &mut Vec<f64> {
        match spin {
            Spin::Medium => &mut self.medium,
            Spin::Solid => &mut self.solid,
            Spin::Some(_) => &mut self.cell
        }
    }
}

fn clat_alat(
    cell_lattice_path: &str,
    act_lattice_path: &str,
    width: usize,
    height: usize
) -> anyhow::Result<(Lattice<Spin>, Lattice<u32>)> {
    let rect = Rect::new((0, 0).into(), (width, height).into());
    let clat = IoManager::read_cell_lattice(
        cell_lattice_path,
        rect.clone()
    ).with_context(|| "while reading cell lattice")?;
    let alat = IoManager::read_lattice_u32(
        act_lattice_path,
        rect.clone()
    ).with_context(|| "while reading act lattice")?;
    Ok((clat, alat))
}