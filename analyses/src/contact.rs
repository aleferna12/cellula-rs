use std::collections::{HashMap, HashSet};
use anyhow::{bail, Context};
use cellulars_lib::constants::CellIndex;
use cellulars_lib::lattice::Lattice;
use cellulars_lib::positional::boundaries::{Boundary, FixedBoundary};
use cellulars_lib::positional::neighbourhood::{MooreNeighbourhood, Neighbourhood};
use cellulars_lib::positional::pos::Pos;
use cellulars_lib::positional::rect::Rect;
use cellulars_lib::spin::Spin;
use num::NumCast;
use polars::polars_utils::float::IsFloat;
use polars::prelude::AnyValue;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;

#[pymodule]
pub mod contact {
    #[pymodule_export]
    use super::geom_act;
    #[pymodule_export]
    use super::kernel_act;
    #[pymodule_export]
    use super::local_act;
    #[pymodule_export]
    use super::neighbour_map;
}

#[pyfunction]
pub fn local_act(
    clat: PyDataFrame,
    alat: PyDataFrame
) -> PyResult<ActVecs> {
    let clat = into_spin_lat(clat).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let alat = into_lat::<u32>(alat).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let neighborhood = MooreNeighbourhood::new(1);
    let bound = FixedBoundary::new(Rect::new(
        clat.rect.min.to_isize(),
        clat.rect.max.to_isize()
    ));
    let mut act_vecs = ActVecs::default();
    for pos in clat.iter_positions() {
        let spin = clat[pos];
        let Spin::Some(_) = spin else {
            continue;
        };
        for neigh_spin in filter_neighs(pos, &clat, &neighborhood, &bound).unwrap() {
            let vec = act_vecs.vec_ref_mut(neigh_spin);
            vec.push(alat[pos] as f64);
        }
    }

    Ok(act_vecs)
}

#[pyfunction]
pub fn geom_act(
    clat: PyDataFrame,
    alat: PyDataFrame
) -> PyResult<ActVecs> {
    let clat = into_spin_lat(clat).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let alat = into_lat::<u32>(alat).map_err(|e| PyValueError::new_err(e.to_string()))?;
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
            let act_diff = neigh_act - act;
            vec.push(act_diff);
        }
    }

    Ok(act_vecs)
}

#[pyfunction]
pub fn kernel_act(
    clat: PyDataFrame,
    alat: PyDataFrame,
    radius: u8,
    geom: bool
) -> PyResult<(ActVecs, Vec<f64>)> {
    let clat = into_spin_lat(clat).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let alat = into_lat::<u32>(alat).map_err(|e| PyValueError::new_err(e.to_string()))?;

    let mut klat = Lattice::new(clat.rect.clone());
    let neighborhood = MooreNeighbourhood::new(1);
    let kernel = MooreNeighbourhood::new(radius);
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

        let mut act_sum = alat[pos];
        let mut neigh_count = 1;
        for neigh in kernel.neighbours(pos.to_isize()) {
            let Some(valid_neigh) = bound.valid_pos(neigh) else {
                continue;
            };
            let lat_neigh = valid_neigh.to_usize();
            let neigh_spin = clat[lat_neigh];
            if let Spin::Some(neighbor_index) = neigh_spin && cell_index == neighbor_index {
                let neigh_act = alat[lat_neigh];
                if geom {
                    act_sum *= neigh_act;
                } else {
                    act_sum += neigh_act;
                }
                neigh_count += 1;
            }
        }

        let act = if geom {
            (act_sum as f64).powf( 1. / neigh_count as f64)
        } else {
            act_sum as f64 / neigh_count as f64
        };
        klat[pos] = act;

        for neigh_spin in filter_neighs(pos, &clat, &neighborhood, &bound).unwrap() {
            let vec = act_vecs.vec_ref_mut(neigh_spin);
            vec.push(klat[pos]);
        }
    }
    Ok((act_vecs, klat.as_array().into()))
}

#[pyfunction]
pub fn neighbour_map(
    clat: PyDataFrame,
    include_self: bool
) -> PyResult<HashMap<String, HashSet<String>>> {
    let clat = into_spin_lat(clat).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let mut neigh_map = HashMap::new();
    let neighborhood = MooreNeighbourhood::new(1);
    let bound = FixedBoundary::new(Rect::new(
        clat.rect.min.to_isize(),
        clat.rect.max.to_isize()
    ));
    for pos in clat.iter_positions() {
        let spin = clat[pos];
        let entry = neigh_map.entry(spin_to_str(spin)).or_insert_with(HashSet::new);
        for neigh in neighborhood.neighbours(pos.to_isize()) {
            let Some(valid_neigh) = bound.valid_pos(neigh) else {
                continue;
            };
            let lat_neigh = valid_neigh.to_usize();
            let neigh_spin = clat[lat_neigh];
            if !include_self && spin == neigh_spin {
                continue;
            }
            entry.insert(spin_to_str(neigh_spin));
        }
    }

    Ok(neigh_map)
}

fn filter_neighs(
    pos: Pos<usize>,
    clat: &Lattice<Spin>,
    neighborhood: &impl Neighbourhood,
    bound: &impl Boundary<Coord = isize>
) -> Option<impl Iterator<Item = Spin>> {
    let Spin::Some(cell_index) = clat[pos] else {
        return None;
    };
    Some(neighborhood.neighbours(pos.to_isize()).filter_map(move |neigh| {
        let Some(valid_neigh) = bound.valid_pos(neigh) else {
            return None;
        };
        let neigh_spin = clat[valid_neigh.to_usize()];
        if let Spin::Some(neighbor_index) = neigh_spin && cell_index == neighbor_index {
            return None;
        }
        Some(neigh_spin)
    }))
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

fn into_spin_lat(df: PyDataFrame) -> anyhow::Result<Lattice<Spin>> {
    let df = df.0;
    let rect = Rect::new((0, 0).into(), (df.width(), df.height()).into());
    let mut lat = Lattice::new(rect);
    for pos in lat.iter_positions() {
        let col = &df[pos.x];
        let maybe_val = col.get(pos.y)?;
        match maybe_val {
            AnyValue::String(val) => {
                let spin = str_to_spin(val)?;
                lat[pos] = spin;
            },
            AnyValue::Null => bail!("cell lattice contains null values"),
            _ => bail!("cell lattice contains invalid value")
        }
    }
    Ok(lat)
}

fn into_lat<T: Clone + Default + NumCast + IsFloat>(df: PyDataFrame) -> anyhow::Result<Lattice<T>> {
    let df = df.0;
    let rect = Rect::new((0, 0).into(), (df.width(), df.height()).into());
    let mut lat = Lattice::new(rect);
    for pos in lat.iter_positions() {
        let x = df[pos.x].get(pos.y)?.try_extract::<T>()?;
        lat[pos] = x;
    }
    Ok(lat)
}

fn spin_to_str(spin: Spin) -> String {
    match spin {
        Spin::Solid => String::from("s"),
        Spin::Medium => String::from("m"),
        Spin::Some(cell_index) => cell_index.to_string(),
    }
}

fn str_to_spin(s: &str) -> anyhow::Result<Spin> {
    Ok(match s {
        "s" => Spin::Solid,
        "m" => Spin::Medium,
        _ => {
            let cell_index = s.parse::<CellIndex>().with_context(|| {
                format!("lattice contains invalid value {s}")
            })?;
            Spin::Some(cell_index)
        },
    })
}