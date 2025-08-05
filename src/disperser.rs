use std::collections::HashMap;
use std::ptr;
use by_address::ByAddress;
use crate::constants::Spin;
use crate::pond::Pond;
use crate::selector::Selector;
use rustworkx_core::connectivity::connected_components;

pub trait Disperser {
    fn disperse<'p>(&mut self, dispersable: &'p [Pond]) -> Vec<DispersionEvent<'p>>;
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
    fn disperse<'p>(&mut self, dispersable: &'p [Pond]) -> Vec<DispersionEvent<'p>> {
        let selected = self.selector.select(dispersable);
        let mut selected_count = HashMap::<
            ByAddress<&Pond>,
            Vec<ByAddress<&Pond>>
        >::new();
        for d in dispersable {
            let descendents = selected.iter()
                .filter_map(|s| {
                    if ptr::eq(d, *s) {
                        Some(ByAddress(*s))
                    } else {
                        None
                    }
                })
                .collect();
            selected_count.insert(ByAddress(d), descendents);
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
                let prop_spins: Vec<_> = subgraphs[prop_index.0]
                    .iter()
                    .map(|s| { s.index() as Spin })
                    .collect();
                for desc in descendents {
                    events.push(DispersionEvent {
                        from: &pond,
                        to: &desc,
                        spins: prop_spins.clone()
                    })
                }
            }
        }
        events
    }
}