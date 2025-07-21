use crate::genome::CellType::{Divide, Migrate};
use crate::genome::GrnError::{MissingNode, WrongNodeType};
use crate::genome::GrnGene::{Input, Output, Regulatory};
use petgraph::graph::{DiGraph, NodeIndex};
use std::ops::Range;
use thiserror::Error;

pub trait Genome {
    fn attempt_mutate(&mut self) -> bool;
}

pub struct Grn {
    graph: DiGraph<GrnGene, f32>,
    input_ids: Range<u32>,
    regulatory_ids: Range<u32>,
    output_ids: Range<u32>
}

impl Grn {
    pub fn new(input_scales: &[f32], n_regulatory: u32, n_outputs: u32) -> Self {
        let n_inputs = input_scales.len() as u32;
        let mut grn = Self {
            graph: DiGraph::new(),
            input_ids: 0..n_inputs,
            regulatory_ids: (n_inputs + n_outputs)..(n_inputs + n_outputs + n_regulatory),
            output_ids: n_inputs..(n_inputs + n_outputs),
        };
        for input in grn.input_ids.clone() {
            grn.graph.add_node(Input {
                value: 0.,
                scale:input_scales[input as usize]
            });
        }
        for _ in grn.regulatory_ids.clone() {
            grn.graph.add_node(Regulatory {
                threshold: 0.,
                active: false,
                activating: false
            });
        }
        for _ in grn.output_ids.clone() {
            grn.graph.add_node(Output {
                threshold: 0.,
                active: false
            });
        }

        for reg in grn.regulatory_ids.clone() {
            for input in grn.input_ids.clone() {
                grn.graph.add_edge(input.into(), reg.into(), 0.);
            }
            for reg2 in grn.regulatory_ids.clone() {
                grn.graph.add_edge(reg.into(), reg2.into(), 0.);
            }
            for output in grn.output_ids.clone() {
                grn.graph.add_edge(reg.into(), output.into(), 0.);
            }
        }
        grn
    }

    pub fn set_input_value(&mut self, index: NodeIndex, new_value: f32) -> Result<(), GrnError> {
        let node = self.graph.node_weight_mut(index).ok_or(MissingNode)?;
        match node {
            Input {value, scale: _ } => {
                *value = new_value;
                Ok(())
            },
            _ => Err(WrongNodeType)
        }
    }
    
    pub fn get_output_state(&self, index: NodeIndex) -> Result<bool, GrnError> {
        let node = self.graph.node_weight(index).ok_or(MissingNode)?;
        match node {
            Output { threshold: _, active } => {
                Ok(*active)
            },
            _ => Err(WrongNodeType)
        }
    }
}

#[derive(Error, Debug)]
pub enum GrnError {
    #[error("The requested node does not exist in the GRN")]
    MissingNode,
    #[error("The requested node is not the right type (perhaps you used the wrong index?)")]
    WrongNodeType
}

pub enum GrnGene {
    Input {
        value: f32,
        scale: f32
    },
    Regulatory {
        threshold: f32,
        active: bool,
        activating: bool
    },
    Output {
        threshold: f32,
        active: bool
    }
}

pub struct SpecializedGrn {
    grn: Grn
}

impl SpecializedGrn {
    pub fn new(light_scale: f32, size_scale: f32) -> Self {
        Self {
            grn: Grn::new(&[light_scale, size_scale], 2, 1)
        }
    }
    
    pub fn set_light(&mut self, value: u32) {
        if let Err(e) = self.grn.set_input_value(0.into(), value as f32) {
            panic!("Failed to set light with error `{e}`");
        }
    }
    
    pub fn set_size(&mut self, value: u32) {
        if let Err(e) = self.grn.set_input_value(1.into(), value as f32) {
            panic!("Failed to set size with error `{e}`");
        }
    }
    
    pub fn get_cell_type(&self) -> CellType {
        match self.grn.get_output_state(2.into()) { 
            Ok(state) => match state {
                false => Migrate,
                true => Divide
            },
            Err(e) => panic!("Failed to get cell type with error `{e}`")
        }
    }
}

pub enum CellType {
    Migrate,
    Divide
}