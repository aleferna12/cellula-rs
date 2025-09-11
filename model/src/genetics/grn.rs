use std::collections::HashSet;
use crate::genetics::genome::Genome;
use crate::genetics::grn::GrnGeneType::{Input, Output, Regulatory};
use rand::Rng;
use rand_distr::Distribution;
use rand_distr::Normal;
use rustworkx_core::petgraph::prelude::*;
use rustworkx_core::petgraph::visit::IntoNodeReferences;
use serde::{Deserialize, Serialize};
use thiserror::Error;

// TODO: make it so Grn can take any Distribution
#[derive(Clone, Debug)]
pub struct Grn<const I: usize, const O: usize> {
    graph: DiGraph<GrnGeneType, f32>,
    input_ids: [NodeIndex; I],
    output_ids: [NodeIndex; O],
    regulatory_ids: Box<[NodeIndex]>,
    pub mut_rate: f32,
    pub mut_distr: Normal<f32>,
    pub input_signals: [f32; I]
}

impl<const I: usize, const O: usize> Grn<I, O> {
    pub fn from_graph(
        graph: DiGraph<GrnGeneType, f32>,
        mut_rate: f32,
        mut_std: f32
    ) -> Result<Self, GrnError> {
        let mut maybe_inputs = [None; I];
        let mut input_signals = [0.0; I];
        let mut input_counter = 0;

        let mut maybe_outputs = [None; O];
        let mut output_counter = 0;

        let mut regulatory_ids = vec![];
        for (nx, node) in graph.node_references() {
            match node {
                Input(gene) => {
                    if input_counter >= I {
                        return Err(GrnError::TooManyInputs);
                    }
                    maybe_inputs[input_counter] = Some(nx);
                    input_signals[input_counter] = gene.signal;
                    input_counter += 1;
                }
                Regulatory(_) => {
                    regulatory_ids.push(nx);
                }
                Output(_) => {
                    if output_counter >= O {
                        return Err(GrnError::TooManyOutputs);
                    }
                    maybe_outputs[output_counter] = Some(nx);
                    output_counter += 1;
                }
            }
        }

        if input_counter < I {
            return Err(GrnError::TooFewInputs);
        }
        if output_counter < O {
            return Err(GrnError::TooFewOutputs);
        }

        fn transform_ids<const S: usize>(
            maybe: [Option<NodeIndex>; S],
            id_set: &mut HashSet<NodeIndex>
        ) -> Result<[NodeIndex; S], GrnError> {
            let mut transformed = [NodeIndex::default(); S];
            for (i, maybe_id) in maybe.iter().enumerate() {
                let id = maybe_id.unwrap();
                if !id_set.insert(id) {
                    return Err(GrnError::RepeatedId(id.index()));
                }
                transformed[i] = id;
            }
            Ok(transformed)
        }

        let mut id_set = HashSet::new();
        let input_ids = transform_ids(maybe_inputs, &mut id_set)?;
        let output_ids = transform_ids(maybe_outputs, &mut id_set)?;

        Ok(Self {
            graph,
            mut_rate,
            input_ids,
            output_ids,
            input_signals,
            regulatory_ids: regulatory_ids.into_boxed_slice(),
            mut_distr: Normal::new(0., mut_std).expect("invalid `mut_std`")
        })
    }

    pub fn from_sampler(
        input_scales: [f32; I],
        n_regulatory: usize,
        mut_rate: f32,
        mut_std: f32,
        mut sampler: impl FnMut() -> f32
    ) -> Self {
        let mut grn = Grn {
            graph: DiGraph::new(),
            input_ids: core::array::from_fn(NodeIndex::new),
            output_ids: core::array::from_fn(|i| NodeIndex::new(i + I)),
            regulatory_ids: ((I + O)..(I + O + n_regulatory)).map(NodeIndex::new).collect(),
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
                grn.graph.add_edge(input, reg, sampler());
            }
            for reg2 in grn.regulatory_ids.iter().copied() {
                grn.graph.add_edge(reg, reg2, sampler());
            }
            for output in grn.output_ids.iter().copied() {
                grn.graph.add_edge(reg, output, sampler());
            }
        }
        grn
    }
    
    pub fn empty() -> Self {
        Self::from_sampler(
            [0.; I],
            0,
            0.,
            0.,
            || 0.
        )
    }
    
    pub fn graph(&self) -> &DiGraph<GrnGeneType, f32> {
        &self.graph
    }

    pub fn nth_input_gene(&self, index: usize) -> &InputGene {
        self.get_input_gene(self.input_ids[index])
    }

    pub fn nth_regulatory_gene(&self, index: usize) -> &RegulatoryGene {
        self.get_regulatory_gene(self.regulatory_ids[index])
    }

    pub fn nth_output_gene(&self, index: usize) -> &OutputGene {
        self.get_output_gene(self.output_ids[index])
    }

    fn get_input_gene(&self, index: NodeIndex) -> &InputGene {
        match self.graph.node_weight(index) {
            None => panic!("node index does not exist"),
            Some(grn_gene) => match grn_gene {
                Input(gene) => gene,
                _ => panic!("node index is wrong gene type (expected `Input`, received `{grn_gene:?}`)")
            }
        }
    }

    fn get_input_gene_mut(&mut self, index: NodeIndex) -> &mut InputGene {
        match self.graph.node_weight_mut(index) {
            None => panic!("node index does not exist"),
            Some(grn_gene) => match grn_gene {
                Input(gene) => gene,
                _ => panic!("node index is wrong gene type (expected `Input`, received `{grn_gene:?}`)")
            }
        }
    }

    fn get_regulatory_gene(&self, index: NodeIndex) -> &RegulatoryGene {
        match self.graph.node_weight(index) {
            None => panic!("node index does not exist"),
            Some(grn_gene) => match grn_gene {
                Regulatory(gene) => gene,
                _ => panic!("node index is wrong gene type (expected `Regulatory`, received `{grn_gene:?}`)")
            }
        }
    }

    fn get_regulatory_gene_mut(&mut self, index: NodeIndex) -> &mut RegulatoryGene {
        match self.graph.node_weight_mut(index) {
            None => panic!("node index does not exist"),
            Some(grn_gene) => match grn_gene {
                Regulatory(gene) => gene,
                _ => panic!("node index is wrong gene type (expected `Regulatory`, received `{grn_gene:?}`)")
            }
        }
    }

    fn get_output_gene(&self, index: NodeIndex) -> &OutputGene {
        match self.graph.node_weight(index) {
            None => panic!("node index does not exist"),
            Some(grn_gene) => match grn_gene {
                Output(gene) => gene,
                _ => panic!("node index is wrong gene type (expected `Output`, received `{grn_gene:?}`)")
            }
        }
    }

    fn get_output_gene_mut(&mut self, index: NodeIndex) -> &mut OutputGene {
        match self.graph.node_weight_mut(index) {
            None => panic!("node index does not exist"),
            Some(grn_gene) => match grn_gene {
                Output(gene) => gene,
                _ => panic!("node index is wrong gene type (expected `Output`, received `{grn_gene:?}`)")
            }
        }
    }

    fn compute_activation_from_inputs(&self, reg: NodeIndex) -> f32 {
        self.input_ids
            .iter()
            .copied()
            .map(|inp_ind| {
                let input_gene = self.get_input_gene(inp_ind);
                let edge = self.graph.find_edge(inp_ind, reg).expect("missing edge");
                input_gene.signal * self.graph.edge_weight(edge).copied().unwrap()
            })
            .sum()
    }

    fn compute_activation_from_regulators(&self, target: NodeIndex) -> f32 {
        self.regulatory_ids
            .iter()
            .copied()
            .map(|reg_ind| {
                let reg_gene = self.get_regulatory_gene(reg_ind);
                if reg_gene.active {
                    let edge = self.graph.find_edge(reg_ind, target).expect("missing edge");
                    self.graph.edge_weight(edge).copied().unwrap()
                } else {
                    0.
                }
            })
            .sum()
    }
}

impl<const I: usize, const O: usize> Genome for Grn<I, O> {
    fn attempt_mutate(&mut self, rng: &mut impl Rng) -> bool {
        let mutated = false;
        for reg in self.regulatory_ids.clone() {
            if rng.random_bool(self.mut_rate as f64) {
                self.get_regulatory_gene_mut(reg).threshold += self.mut_distr.sample(rng);
            }

            for from in self.input_ids.iter().copied().chain(self.regulatory_ids.iter().copied()) {
                let edge_index = self.graph.find_edge(
                    from,
                    reg
                ).expect("missing edge");
                let edge = self.graph.edge_weight_mut(edge_index).unwrap();
                if rng.random_bool(self.mut_rate as f64) {
                    *edge += self.mut_distr.sample(rng);
                }
            }

            for out in self.output_ids.iter().copied() {
                let edge_index = self.graph.find_edge(
                    reg,
                    out
                ).expect("missing edge");
                let edge = self.graph.edge_weight_mut(edge_index).unwrap();
                *edge += self.mut_distr.sample(rng);
            }
        }

        for out in self.output_ids {
            if rng.random_bool(self.mut_rate as f64) {
                self.get_output_gene_mut(out).threshold += self.mut_distr.sample(rng);
            }
        }
        mutated
    }

    fn update_expression(&mut self) {
        for (i, inp) in self.input_ids.into_iter().enumerate() {
            let signal = self.input_signals[i];
            let inp_gene = self.get_input_gene_mut(inp);
            inp_gene.signal = signal;
        }

        for reg in self.regulatory_ids.clone() {
            let act_influx = self.compute_activation_from_inputs(reg)
                + self.compute_activation_from_regulators(reg);
            let reg_gene = self.get_regulatory_gene_mut(reg);
            reg_gene.activating = act_influx > reg_gene.threshold;
        }
        for reg in self.regulatory_ids.clone() {
            let reg_gene = self.get_regulatory_gene_mut(reg);
            reg_gene.active = reg_gene.activating;
        }

        for out in self.output_ids {
            let act_influx = self.compute_activation_from_regulators(out);
            let out_gene = self.get_output_gene_mut(out);
            out_gene.active = act_influx > out_gene.threshold;
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct InputGene {
    pub signal: f32,
    pub scale: f32
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RegulatoryGene {
    pub threshold: f32,
    pub active: bool,
    activating: bool
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct OutputGene {
    pub threshold: f32,
    pub active: bool
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "kebab-case")]
pub enum GrnGeneType {
    Input(InputGene),
    Regulatory(RegulatoryGene),
    Output(OutputGene)
}

#[derive(Serialize, Deserialize)]
pub struct EdgeWeight {
    pub weight: f32,
}

#[derive(Error, Debug, Clone)]
pub enum GrnError {
    #[error("too many input genes in the graph")]
    TooManyInputs,
    #[error("too few input genes in the graph")]
    TooFewInputs,
    #[error("too many output genes in the graph")]
    TooManyOutputs,
    #[error("too few output genes in the graph")]
    TooFewOutputs,
    #[error("gene with id `{0}` appeared multiple times")]
    RepeatedId(usize)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256StarStar;

    fn fixed_sampler() -> impl FnMut() -> f32 {
        let mut val = 0.0;
        move || {
            val += 1.0;
            val
        }
    }

    // ---------- GRN Tests ----------

    #[test]
    fn graph_has_expected_nodes_and_edges() {
        let grn: Grn<2, 1> = Grn::from_sampler([1.0, 2.0], 1, 0.1, 0.5, fixed_sampler());

        // 2 input + 1 output + 1 regulatory
        assert_eq!(grn.graph.node_count(), 4);

        // Each regulatory node should have edges from inputs, from itself, and to outputs
        let reg = grn.regulatory_ids[0];
        assert!(grn.graph.find_edge(grn.input_ids[0], reg).is_some());
        assert!(grn.graph.find_edge(grn.input_ids[1], reg).is_some());
        assert!(grn.graph.find_edge(reg, reg).is_some());
        assert!(grn.graph.find_edge(reg, grn.output_ids[0]).is_some());
    }

    #[test]
    fn nth_gene_accessors() {
        let grn: Grn<2, 1> = Grn::from_sampler([1.5, 2.5], 1, 0.1, 0.5, fixed_sampler());

        assert_eq!(grn.nth_input_gene(0).scale, 1.5);
        assert_eq!(grn.nth_input_gene(1).scale, 2.5);
        assert!(grn.nth_output_gene(0).threshold > 0.0);
        assert!(grn.nth_regulatory_gene(0).threshold > 0.0);
    }

    #[test]
    fn update_expression_sets_signals_and_activation() {
        let mut grn: Grn<1, 1> = Grn::from_sampler([1.0], 1, 0.1, 0.5, || 1.0);

        // Give input a signal so activation > threshold
        grn.input_signals[0] = 10.0;
        grn.update_expression();

        let reg_gene = grn.nth_regulatory_gene(0);
        assert!(reg_gene.active);
    }

    #[test]
    fn attempt_mutate_changes_thresholds_and_edges() {
        let mut grn: Grn<1, 1> = Grn::from_sampler([1.0], 1, 1.0, 0.5, || 0.0);
        let mut rng = Xoshiro256StarStar::seed_from_u64(42);

        // Snapshot values before
        let reg_threshold_before = grn.nth_regulatory_gene(0).threshold;
        let out_threshold_before = grn.nth_output_gene(0).threshold;

        // Collect all edge weights before mutation
        let edges_before: Vec<_> = grn.graph.edge_weights().copied().collect();

        grn.attempt_mutate(&mut rng);

        // Thresholds should change
        assert_ne!(grn.nth_regulatory_gene(0).threshold, reg_threshold_before);
        assert_ne!(grn.nth_output_gene(0).threshold, out_threshold_before);

        // At least one edge weight should have changed
        let edges_after: Vec<_> = grn.graph.edge_weights().copied().collect();
        assert!(edges_before.iter().zip(edges_after.iter()).all(|(a, b)| a != b));
    }
}
