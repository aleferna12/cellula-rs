use crate::constants::Spin;
use crate::ecology::selector::PreservesOrder;
use crate::pond::Pond;
use rustworkx_core::connectivity::connected_components;

pub trait Disperser {
    fn disperse(&mut self, dispersable: &[Pond]) -> Vec<DispersionEvent>;
}

// This can become a trait in the future if needed
pub struct DispersionEvent {
    from: usize,
    to: usize,
    spins: Vec<Spin>,
}

pub struct SelectiveDispersion<S> {
    selector: S
}

impl<S> SelectiveDispersion<S> {
    pub fn get_prop_spins(pond: &Pond) -> Vec<Spin> {
        let neighs_graph = pond.env.build_neighbours_graph();
        let subgraphs = connected_components(&neighs_graph);
        // There is only one cluster
        if subgraphs.len() <= 1 {
            return vec![];
        }

        let prop_index = subgraphs
            .iter()
            .map(|s| s.len())
            .enumerate()
            .min_by(|a, b| a.1.cmp(&b.1))
            .unwrap();

        subgraphs[prop_index.0]
            .iter()
            .map(|s| { s.index() as Spin })
            .collect()
    }
}

impl<S: PreservesOrder> Disperser for SelectiveDispersion<S> {
    fn disperse(&mut self, dispersable: &[Pond]) -> Vec<DispersionEvent> {
        let selected = self.selector.select(dispersable);
        let mut events = vec![];
        // Lazily evaluated
        let mut props = vec![None; selected.len()];
        for (i, parent) in selected.into_iter().enumerate() {
            if i == parent {
                continue;
            }

            let prop = match props[parent].as_ref() {
                Some(prop) => prop,
                None => {
                    props[parent] = Self::get_prop_spins(&dispersable[parent]).into();
                    props[parent].as_ref().unwrap()
                }
            };

            if !prop.is_empty() {
                events.push(DispersionEvent {
                    from: parent,
                    to: i,
                    spins: prop.clone(),
                })
            }
        }
        events
    }
}