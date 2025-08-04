use crate::genome::GrnError::{MissingNode, WrongNodeType};
use crate::genome::GrnGeneType::{Input, Output, Regulatory};
use rand::Rng;
use rand_distr::Distribution;
use rand_distr::Normal;
use rustworkx_core::petgraph::prelude::{DiGraph, NodeIndex};
use thiserror::Error;

pub trait Genome {
    fn attempt_mutate(&mut self, rng: &mut impl Rng) -> bool;
    fn update_expression(&mut self);
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

    fn update_expression(&mut self) {
        self.counter += 1;
        if self.counter > self.period_updates {
            match self.cell_type {
                CellType::Migrate => self.cell_type = CellType::Divide,
                CellType::Divide => self.cell_type = CellType::Migrate
            }
            self.counter = 0;
        }
    }
}

/// Provides a generic, versatile GRN implementation that is specialised in `Grn`.
#[derive(Clone, Debug)]
pub struct Grn<const I: usize, const O: usize> {
    graph: DiGraph<GrnGeneType, f32>,
    input_ids: [usize; I],
    output_ids: [usize; I],
    regulatory_ids: Vec<usize>,
    mut_rate: f32,
    mut_distr: Normal<f32>,
    pub input_signals: [f32; I]
}

impl<const I: usize, const O: usize> Grn<I, O> {
    pub fn new(
        input_scales: [f32; I],
        n_regulatory: usize,
        mut_rate: f32,
        mut_std: f32,
        mut sampler: impl FnMut() -> f32
    ) -> Self {
        let mut grn = Grn {
            graph: DiGraph::new(),
            input_ids: core::array::from_fn(|i| i),
            output_ids: core::array::from_fn(|i| i + I),
            regulatory_ids: ((I + O)..(I + O + n_regulatory)).collect(),
            input_signals: [0.; I],
            mut_rate,
            mut_distr: Normal::new(0., mut_std).expect("invalid `mut_std`")
        };

        for scale in input_scales {
            grn.graph.add_node(Input(InputGene { signal: 0., scale }));
        }
        for _ in 0..O {
            grn.graph.add_node(Output(OutputGene {
                threshold: sampler(),
                active: false
            }));
        }
        for _ in 0..n_regulatory {
            grn.graph.add_node(Regulatory(RegulatoryGene {
                threshold: sampler(),
                active: false,
                activating: false
            }));
        }

        for reg in grn.regulatory_ids.iter().copied() {
            for input in grn.input_ids.iter().copied() {
                grn.graph.add_edge(NodeIndex::new(input), NodeIndex::new(reg), sampler());
            }
            for reg2 in grn.regulatory_ids.iter().copied() {
                grn.graph.add_edge(NodeIndex::new(reg), NodeIndex::new(reg2), sampler());
            }
            for output in grn.output_ids.iter().copied() {
                grn.graph.add_edge(NodeIndex::new(reg), NodeIndex::new(output), sampler());
            }
        }
        grn
    }

    pub fn get_cell_type(&self) -> CellType {
        match self
            .get_output_gene(NodeIndex::new(self.output_ids[0]))
            .expect("invalid GRN architecture").active {
            false => CellType::Migrate,
            true => CellType::Divide
        }
    }

    fn get_input_gene(&self, index: NodeIndex) -> Result<&InputGene, GrnError> {
        let node = self.graph.node_weight(index).ok_or(MissingNode)?;
        match node {
            Input(gene) => {
                Ok(gene)
            },
            _ => Err(WrongNodeType)
        }
    }

    fn get_input_gene_mut(&mut self, index: NodeIndex) -> Result<&mut InputGene, GrnError> {
        let node = self.graph.node_weight_mut(index).ok_or(MissingNode)?;
        match node {
            Input(gene) => {
                Ok(gene)
            },
            _ => Err(WrongNodeType)
        }
    }

    fn get_regulatory_gene(&self, index: NodeIndex) -> Result<&RegulatoryGene, GrnError> {
        let node = self.graph.node_weight(index).ok_or(MissingNode)?;
        match node {
            Regulatory(gene) => {
                Ok(gene)
            },
            _ => Err(WrongNodeType)
        }
    }

    fn get_regulatory_gene_mut(&mut self, index: NodeIndex) -> Result<&mut RegulatoryGene, GrnError> {
        let node = self.graph.node_weight_mut(index).ok_or(MissingNode)?;
        match node {
            Regulatory(gene) => {
                Ok(gene)
            },
            _ => Err(WrongNodeType)
        }
    }

    fn get_output_gene(&self, index: NodeIndex) -> Result<&OutputGene, GrnError> {
        let node = self.graph.node_weight(index).ok_or(MissingNode)?;
        match node {
            Output(gene) => {
                Ok(gene)
            },
            _ => Err(WrongNodeType)
        }
    }

    fn get_output_gene_mut(&mut self, index: NodeIndex) -> Result<&mut OutputGene, GrnError> {
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
            let inp_ind = NodeIndex::new(inp);
            let input_gene = self.get_input_gene(inp_ind).expect("missing input gene");
            let edge = self.graph.find_edge(inp_ind, reg)?;
            act += input_gene.signal * self.graph.edge_weight(edge).copied()?
        }
        Some(act)
    }

    fn compute_activation_from_regulators(&self, target: NodeIndex) -> Option<f32> {
        let mut act = 0.;
        for reg in self.regulatory_ids.iter().copied() {
            let reg_ind = NodeIndex::new(reg);
            let reg_gene = self.get_regulatory_gene(reg_ind).expect("missing regulatory gene");
            act += if reg_gene.active {
                let edge = self.graph.find_edge(reg_ind, target)?;
                self.graph.edge_weight(edge).copied()?
            } else {
                0.
            }
        }
        Some(act)
    }
}

impl<const I: usize, const O: usize> Genome for Grn<I, O> {
    fn attempt_mutate(&mut self, rng: &mut impl Rng) -> bool {
        let mutated = false;
        let reg_miss = "Missing regulatory gene";
        for reg in self.regulatory_ids.clone() {
            if rng.random_bool(self.mut_rate as f64) {
                self.get_regulatory_gene_mut(NodeIndex::new(reg))
                    .expect(reg_miss)
                    .threshold += self.mut_distr.sample(rng);
            }

            for from in self.input_ids.iter().copied().chain(self.regulatory_ids.iter().copied()) {
                let edge_index = self.graph.find_edge(
                    NodeIndex::new(from),
                    NodeIndex::new(reg)
                ).expect("missing edge");
                let edge = self.graph.edge_weight_mut(edge_index).unwrap();
                if rng.random_bool(self.mut_rate as f64) {
                    *edge += self.mut_distr.sample(rng);
                }
            }

            for out in self.output_ids.iter().copied() {
                let edge_index = self.graph.find_edge(
                    NodeIndex::new(reg),
                    NodeIndex::new(out)
                ).expect("missing edge");
                let edge = self.graph.edge_weight_mut(edge_index).unwrap();
                *edge += self.mut_distr.sample(rng);
            }
        }

        for out in self.output_ids {
            if rng.random_bool(self.mut_rate as f64) {
                self.get_output_gene_mut(NodeIndex::new(out))
                    .expect(reg_miss)
                    .threshold += self.mut_distr.sample(rng);
            }
        }
        mutated
    }

    fn update_expression(&mut self) {
        for (i, inp) in self.input_ids.into_iter().enumerate() {
            let signal = self.input_signals[i];
            let inp_gene = self.get_input_gene_mut(
                NodeIndex::new(inp)
            ).expect("Missing regulatory gene");
            inp_gene.signal = signal;
        }

        let reg_miss = "Missing regulatory gene";
        for reg in self.regulatory_ids.clone() {
            let reg_idx = NodeIndex::new(reg);
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
        for out in self.output_ids {
            let out_ind = NodeIndex::new(out);
            let act_influx = self.compute_activation_from_regulators(out_ind).expect(out_miss);

            let out_gene = self.get_output_gene_mut(out_ind).expect(out_miss);
            out_gene.active = act_influx > out_gene.threshold;
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

#[derive(Clone, Debug)]
pub enum CellType {
    Migrate,
    Divide
}