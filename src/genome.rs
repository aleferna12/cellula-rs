use crate::genome::GrnError::{MissingNode, WrongNodeType};
use crate::genome::GrnGeneType::{Input, Output, Regulatory};
use petgraph::graph::{DiGraph, NodeIndex};
use thiserror::Error;

pub trait Genome {
    fn attempt_mutate(&mut self) -> bool;
    fn update_expression(&mut self);
    fn get_cell_type(&self) -> CellType;
}

#[derive(Clone, Debug)]
pub struct Grn {
    graph: DiGraph<GrnGeneType, f32>,
    input_ids: Vec<u32>,
    regulatory_ids: Vec<u32>,
    output_ids: Vec<u32>
}

// TODO: this might be quite slow, we can implement a StaticGrn that does not support node insertion or deletion
//  for improved performance. The advantage of this implementation is that we can delete or duplicate nodes dynamically
impl Grn {
    pub fn new(input_scales: &[f32], n_regulatory: u32, n_outputs: u32) -> Self {
        let n_inputs = input_scales.len() as u32;
        let mut grn = Grn {
            graph: DiGraph::new(),
            input_ids: Vec::new(),
            regulatory_ids: Vec::new(),
            output_ids: Vec::new(),
        };

        for _ in 0..n_inputs {
            let idx = grn.graph.add_node(Input(InputGene { value: 0. }));
            grn.input_ids.push(idx.index() as u32);
        }
        for _ in 0..n_regulatory {
            let idx = grn.graph.add_node(Regulatory(RegulatoryGene {
                threshold: 0.,
                active: false,
                activating: false
            }));
            grn.regulatory_ids.push(idx.index() as u32);
        }
        for _ in 0..n_outputs {
            let idx = grn.graph.add_node(Output(OutputGene { threshold: 0., active: false }));
            grn.output_ids.push(idx.index() as u32);
        }

        for reg in grn.regulatory_ids.iter().copied() {
            for input in grn.input_ids.iter().copied() {
                grn.graph.add_edge(input.into(), reg.into(), input_scales[input as usize]);
            }
            for reg2 in grn.regulatory_ids.iter().copied() {
                grn.graph.add_edge(reg.into(), reg2.into(), 0.);
            }
            for output in grn.output_ids.iter().copied() {
                grn.graph.add_edge(reg.into(), output.into(), 0.);
            }
        }
        grn
    }

    pub fn get_input_gene(&self, index: NodeIndex) -> Result<&InputGene, GrnError> {
        let node = self.graph.node_weight(index).ok_or(MissingNode)?;
        match node {
            Input(gene) => {
                Ok(gene)
            },
            _ => Err(WrongNodeType)
        }
    }

    pub fn get_input_gene_mut(&mut self, index: NodeIndex) -> Result<&mut InputGene, GrnError> {
        let node = self.graph.node_weight_mut(index).ok_or(MissingNode)?;
        match node {
            Input(gene) => {
                Ok(gene)
            },
            _ => Err(WrongNodeType)
        }
    }

    pub fn get_regulatory_gene(&self, index: NodeIndex) -> Result<&RegulatoryGene, GrnError> {
        let node = self.graph.node_weight(index).ok_or(MissingNode)?;
        match node {
            Regulatory(gene) => {
                Ok(gene)
            },
            _ => Err(WrongNodeType)
        }
    }

    pub fn get_regulatory_gene_mut(&mut self, index: NodeIndex) -> Result<&mut RegulatoryGene, GrnError> {
        let node = self.graph.node_weight_mut(index).ok_or(MissingNode)?;
        match node {
            Regulatory(gene) => {
                Ok(gene)
            },
            _ => Err(WrongNodeType)
        }
    }

    pub fn get_output_gene(&self, index: NodeIndex) -> Result<&OutputGene, GrnError> {
        let node = self.graph.node_weight(index).ok_or(MissingNode)?;
        match node {
            Output(gene) => {
                Ok(gene)
            },
            _ => Err(WrongNodeType)
        }
    }

    pub fn get_output_gene_mut(&mut self, index: NodeIndex) -> Result<&mut OutputGene, GrnError> {
        let node = self.graph.node_weight_mut(index).ok_or(MissingNode)?;
        match node {
            Output(gene) => {
                Ok(gene)
            },
            _ => Err(WrongNodeType)
        }
    }

    fn get_edge_weight(&self, from: NodeIndex, to: NodeIndex) -> Option<f32> {
        let edge = self.graph.find_edge(from, to)?;
        self.graph.edge_weight(edge).copied()
    }

    fn compute_activation_from_inputs(&self, reg: NodeIndex) -> Option<f32> {
        let weights = self.input_ids.iter().copied().map(|inp| {
            let input_gene = self.get_input_gene(inp.into()).expect("Input node missing");
            Some(input_gene.value * self.get_edge_weight(inp.into(), reg)?)
        });
        let mut act = 0.;
        for w in weights {
            act += w?;
        }
        Some(act)
    }

    fn compute_activation_from_regulators(&self, target: NodeIndex) -> Option<f32> {
        let weights = self.regulatory_ids.iter().copied().map(|reg| {
            let reg_gene = self.get_regulatory_gene(reg.into()).expect("Regulatory node missing");
            if reg_gene.active {
                Some(self.get_edge_weight(reg.into(), target)?)
            } else {
                Some(0.)
            }
        });
        let mut act = 0.;
        for w in weights {
            act += w?;
        }
        Some(act)
    }

    fn update_expression(&mut self) {
        let reg_miss = "Regulatory gene missing";
        for reg in self.regulatory_ids.clone() {
            let reg_idx = reg.into();
            let act_influx = self
                .compute_activation_from_inputs(reg_idx)
                .expect(reg_miss)
                + self
                .compute_activation_from_regulators(reg_idx)
                .expect(reg_miss);
            let reg_gene = self.get_regulatory_gene_mut(reg_idx).expect(reg_miss);
            reg_gene.activating = act_influx > reg_gene.threshold;
        }

        let out_miss = "Output gene missing";
        for out in self.output_ids.clone() {
            let out_idx = out.into();
            let act_influx = self.compute_activation_from_regulators(out_idx).expect(out_miss);

            let out_gene = self.get_output_gene_mut(out_idx).expect(out_miss);
            out_gene.active = act_influx > out_gene.threshold;
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

#[derive(Clone)]
#[derive(Debug)]
pub struct InputGene {
    pub value: f32
}

#[derive(Clone)]
#[derive(Debug)]
pub struct RegulatoryGene {
    pub threshold: f32,
    pub active: bool,
    pub activating: bool
}

#[derive(Clone)]
#[derive(Debug)]
pub struct OutputGene {
    pub threshold: f32,
    pub active: bool
}

#[derive(Clone, Debug)]
pub enum GrnGeneType {
    Input(InputGene),
    Regulatory(RegulatoryGene),
    Output(OutputGene)
}

#[derive(Clone)]
#[derive(Debug)]
pub struct SpecialisedGrn {
    grn: Grn
}

impl SpecialisedGrn {
    pub fn new(light_scale: f32, size_scale: f32) -> Self {
        Self {
            grn: Grn::new(&[light_scale, size_scale], 2, 1)
        }
    }
}

impl Genome for SpecialisedGrn {
    fn attempt_mutate(&mut self) -> bool {
        todo!()
    }

    fn update_expression(&mut self) {
        self.grn.update_expression();
    }

    fn get_cell_type(&self) -> CellType {
        let base_index = self.grn.input_ids.len() + self.grn.regulatory_ids.len();
        match self.grn.get_output_gene(NodeIndex::new(base_index)).expect("Missing cell type node").active {
            false => CellType::Migrate,
            true => CellType::Divide
        }
    }
}

#[derive(Clone, Debug)]
pub enum CellType {
    Migrate,
    Divide
}