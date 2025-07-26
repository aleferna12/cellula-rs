use crate::genome::GrnError::{MissingNode, WrongNodeType};
use crate::genome::GrnGeneType::{Input, Output, Regulatory};
use petgraph::graph::{DiGraph, NodeIndex};
use rand::Rng;
use rand_distr::Distribution;
use rand_distr::Normal;
use thiserror::Error;

pub trait Genome {
    fn attempt_mutate(&mut self, rng: &mut impl Rng) -> bool;
    fn update_expression(&mut self, chem_signal: u32);
    fn get_cell_type(&self) -> CellType;
}

/// This is a fake genome that just cycles through cell types.
#[derive(Clone, Debug)]
pub struct MockGenome {
    period_updates: u32,
    counter: u32,
    cell_type: CellType
}

impl MockGenome {
    /// Makes a new `MockGenome` with a specified period.
    ///
    /// `period_updates` is the period for which each cell type will last for.
    /// The unit is the number of `update_expression()` calls, not MCS.
    pub fn new(period_updates: u32) -> Self {
        Self {
            period_updates,
            counter: 0,
            cell_type: CellType::Migrate
        }
    }
}

impl Genome for MockGenome {
    fn attempt_mutate(&mut self, _rng: &mut impl Rng) -> bool {
        false
    }

    fn update_expression(&mut self, _chem_signal: u32) {
        self.counter += 1;
        if self.counter > self.period_updates {
            match self.cell_type {
                CellType::Migrate => self.cell_type = CellType::Divide,
                CellType::Divide => self.cell_type = CellType::Migrate
            }
            self.counter = 0;
        }
    }

    fn get_cell_type(&self) -> CellType {
        self.cell_type.clone()
    }
}

// TODO: this might be quite slow, we can implement a StaticGrn that does not support node insertion or deletion
//  for improved performance. The advantage of this implementation is that we can delete or duplicate nodes dynamically
//  and the object is stored directly as a graph.
/// Provides a generic, versatile GRN implementation that is specialised in `Grn`.
#[derive(Clone, Debug)]
pub struct BaseGrn {
    graph: DiGraph<GrnGeneType, f32>,
    input_ids: Vec<u32>,
    regulatory_ids: Vec<u32>,
    output_ids: Vec<u32>,
    mut_rate: f32,
    mut_distr: Normal<f32>
}

impl BaseGrn {
    pub fn new(
        input_scales: &[f32],
        n_regulatory: u32,
        n_outputs: u32,
        mut_rate: f32,
        mut_std: f32, 
        mut sampler: impl FnMut() -> f32
    ) -> Self {
        let mut grn = BaseGrn {
            graph: DiGraph::new(),
            input_ids: Vec::new(),
            regulatory_ids: Vec::new(),
            output_ids: Vec::new(),
            mut_rate,
            mut_distr: Normal::new(0., mut_std).expect("invalid `mut_std`")
        };

        for scale in input_scales.iter() {
            let idx = grn.graph.add_node(Input(InputGene { signal: 0., scale: *scale }));
            grn.input_ids.push(idx.index() as u32);
        }

        for _ in 0..n_regulatory {
            let idx = grn.graph.add_node(Regulatory(RegulatoryGene {
                threshold: sampler(),
                active: false,
                activating: false
            }));
            grn.regulatory_ids.push(idx.index() as u32);
        }
        for _ in 0..n_outputs {
            let idx = grn.graph.add_node(Output(OutputGene { 
                threshold: sampler(),
                active: false 
            }));
            grn.output_ids.push(idx.index() as u32);
        }

        for reg in grn.regulatory_ids.iter().copied() {
            for input in grn.input_ids.iter().copied() {
                grn.graph.add_edge(input.into(), reg.into(), sampler());
            }
            for reg2 in grn.regulatory_ids.iter().copied() {
                grn.graph.add_edge(reg.into(), reg2.into(), sampler());
            }
            for output in grn.output_ids.iter().copied() {
                grn.graph.add_edge(reg.into(), output.into(), sampler());
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

    fn compute_activation_from_inputs(&self, reg: NodeIndex) -> Option<f32> {
        let mut act = 0.;
        for inp in self.input_ids.iter().copied() {
            let input_gene = self.get_input_gene(inp.into()).expect("missing input gene");
            let edge = self.graph.find_edge(inp.into(), reg)?;
            act += input_gene.signal * self.graph.edge_weight(edge).copied()?
        }
        Some(act)
    }

    fn compute_activation_from_regulators(&self, target: NodeIndex) -> Option<f32> {
        let mut act = 0.;
        for reg in self.regulatory_ids.iter().copied() {
            let reg_gene = self.get_regulatory_gene(reg.into()).expect("missing regulatory gene");
            act += if reg_gene.active {
                let edge = self.graph.find_edge(reg.into(), target)?;
                self.graph.edge_weight(edge).copied()?
            } else {
                0.
            }
        }
        Some(act)
    }

    fn update_expression(&mut self, input_signals: &[f32]) -> Result<(), GrnError> {
        for (i, inp) in self.input_ids.clone().into_iter().enumerate() {
            let inp_gene = self.get_input_gene_mut(inp.into())?;
            inp_gene.signal = input_signals[i];
        }
        
        let reg_miss = "Missing regulatory gene";
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

        let out_miss = "Missing output gene";
        for out in self.output_ids.clone() {
            let out_idx = out.into();
            let act_influx = self.compute_activation_from_regulators(out_idx).expect(out_miss);

            let out_gene = self.get_output_gene_mut(out_idx).expect(out_miss);
            out_gene.active = act_influx > out_gene.threshold;
        }
        Ok(())
    }
    
    fn attempt_mutate(&mut self, rng: &mut impl Rng) -> bool {
        let mutated = false;
        let reg_miss = "Missing regulatory gene";
        for reg in self.regulatory_ids.clone() {
            if rng.random_bool(self.mut_rate as f64) {
                self.get_regulatory_gene_mut(reg.into())
                    .expect(reg_miss)
                    .threshold += self.mut_distr.sample(rng);
            }
            
            for from in self.input_ids.iter().copied().chain(self.regulatory_ids.iter().copied()) {
                let edge_index = self.graph.find_edge(from.into(), reg.into()).expect("missing edge");
                let edge = self.graph.edge_weight_mut(edge_index).unwrap();
                if rng.random_bool(self.mut_rate as f64) {
                    *edge += self.mut_distr.sample(rng);
                }
            }
            
            for out in self.output_ids.iter().copied() {
                let edge_index = self.graph.find_edge(reg.into(), out.into()).expect("missing edge");
                let edge = self.graph.edge_weight_mut(edge_index).unwrap();
                *edge += self.mut_distr.sample(rng);
            }
        }
        mutated
    }
}

/// Specialised `Grn` that handles the correct inputs and outputs.
impl Grn {
    pub fn new(
        chem_scale: f32,
        n_regulatory: u32,
        mut_rate: f32,
        mut_std: f32,
        sampler: impl FnMut() -> f32
    ) -> Self {
        Self {
            grn: BaseGrn::new(
                &[chem_scale],
                n_regulatory, 
                1,
                mut_rate,
                mut_std,
                sampler
            ),
        }
    }
}

impl Genome for Grn {
    fn attempt_mutate(&mut self, rng: &mut impl Rng) -> bool {
        self.grn.attempt_mutate(rng)
    }

    fn update_expression(&mut self, chem_signal: u32) {
        self.grn.update_expression(&[chem_signal as f32]).expect("invalid GRN architecture");
    }

    fn get_cell_type(&self) -> CellType {
        match self.grn
            .get_output_gene(self.grn.output_ids[0].into())
            .expect("invalid GRN architecture").active {
            false => CellType::Migrate,
            true => CellType::Divide
        }
    }
}

#[derive(Error, Debug)]
pub enum GrnError {
    #[error("the requested node does not exist in the GRN")]
    MissingNode,
    #[error("the requested node is not the right type (perhaps you used the wrong index?)")]
    WrongNodeType
}

#[derive(Clone)]
#[derive(Debug)]
pub struct InputGene {
    pub signal: f32,
    pub scale: f32
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
pub struct Grn {
    grn: BaseGrn,
}

#[derive(Clone, Debug)]
pub enum CellType {
    Migrate,
    Divide
}