use std::collections::HashMap;
use std::ptr;
use crate::constants::Spin;
use crate::pond::Pond;
use crate::selector::Selector;
use rustworkx_core::connectivity::connected_components;

pub trait Disperser {
    fn disperse(&mut self, dispersable: &[Pond]) -> Vec<DispersionEvent>;
}

// This can become a trait in the future if needed
pub struct DispersionEvent<'p> {
    from: &'p Pond,
    to: &'p Pond,
    spins: Vec<Spin>,
}

pub struct SelectiveDispersion<S> {
    selector: S
}

impl<S: Selector> Disperser for SelectiveDispersion<S> {
    fn disperse(&mut self, dispersable: &[Pond]) -> Vec<DispersionEvent> {
        let selected = self.selector.select(dispersable);
        let mut selected_count = HashMap::<&Pond, Vec<&Pond>>::new();
        for d in dispersable {
            let descendents = selected.iter()
                .filter(|s| ptr::eq(d, **s))
                .collect();
            selected_count.insert(d, descendents);
        }

        let mut events = vec![];
        for (pond, descendents) in selected_count.into_iter() {
            if descendents.len() < 2 {
                continue
            }

            let neighs_graph = pond.env.build_neighbours_graph();
            let subgraphs = connected_components(&neighs_graph);
            // Is None if the Pond is empty
            if let Some(prop_index) = subgraphs
                .iter()
                .map(|s| { s.len() })
                .enumerate()
                .min_by(|a, b| a.1.cmp(&b.1))  {
                let prop_spins = subgraphs[prop_index.0]
                    .iter()
                    .map(|s| { Spin::from(s) })
                    .collect();
                for desc in descendents {
                    events.push(DispersionEvent {
                        from: pond,
                        to: desc,
                        spins: prop_spins.clone()
                    })
                }
            }
        }
        events
    }
}