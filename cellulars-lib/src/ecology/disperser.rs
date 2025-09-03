use crate::constants::Spin;
use crate::ecology::selector::PreservesOrder;
use crate::pond::Pond;
use rustworkx_core::connectivity::connected_components;

pub trait Disperser {
    fn disperse(&mut self, dispersable: &[Pond]) -> Vec<DispersionEvent>;
}

// This can become a trait in the future if needed
#[derive(Debug)]
pub struct DispersionEvent {
    pub(crate) from: usize,
    pub(crate) to: usize,
    pub(crate) spins: Vec<Spin>,
}

pub struct SelectiveDispersion<S> {
    pub selector: S
}

impl<S> SelectiveDispersion<S> {
    /// Returns at most `n_props` 
    pub fn get_prop_spins(pond: &Pond, n_props: usize) -> Vec<Vec<Spin>> {
        if n_props < 1 {
            return vec![];
        }
        
        let neighs_graph = pond.env.build_neighbours_graph();
        let mut subgraphs = connected_components(&neighs_graph);
        // There is only one cluster
        if subgraphs.len() <= 1 {
            return vec![];
        }

        subgraphs.sort_by(|subgraph1, subgraph2| { 
            subgraph1.len().cmp(&subgraph2.len()) 
        });
        
        subgraphs[0..n_props.min(subgraphs.len() - 1)]
            .iter()
            .map(|subgraph| {
                subgraph.iter().map(|&index| neighs_graph[index]).collect()
            })
            .collect()
    }
}

impl<S: PreservesOrder> Disperser for SelectiveDispersion<S> {
    fn disperse(&mut self, dispersable: &[Pond]) -> Vec<DispersionEvent> {
        let selected = self.selector.select(dispersable);
        let mut prop_counts = vec![0usize; selected.len()];
        for &s in &selected {
            prop_counts[s] += 1;
        }
        let mut props: Vec<_> = prop_counts.into_iter()
            .enumerate()
            .map(|(i, count)| Self::get_prop_spins(&dispersable[i], count.saturating_sub(1)))
            .collect();
        
        selected.into_iter()
            .enumerate()
            .filter_map(|(to, from)| {
                if from == to {
                    None
                } else {
                    props[from].pop().map(|prop| DispersionEvent {
                        from,
                        to,
                        spins: prop
                    })
                }
            })
            .collect()
    }
}