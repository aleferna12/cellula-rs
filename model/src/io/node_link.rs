use crate::genetics::grn::{EdgeWeight, Grn, GrnGeneType};
use anyhow::bail;
use rustworkx_core::petgraph::prelude::{DiGraph, EdgeRef, NodeIndex};
use rustworkx_core::petgraph::visit::IntoNodeReferences;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub struct Node<N> {
    pub id: usize,
    #[serde(flatten)]
    pub data: N,
}

#[derive(Serialize, Deserialize)]
pub struct Link<E> {
    pub source: usize,
    pub target: usize,
    #[serde(flatten)]
    pub data: E,
}

#[derive(Serialize, Deserialize)]
pub struct NodeLinkData<N, E, G> {
    pub directed: bool,
    pub multigraph: bool,
    pub graph: G,
    pub nodes: Vec<Node<N>>,
    pub links: Vec<Link<E>>,
}

impl<const I: usize, const O: usize> From<Grn<I, O>> for NodeLinkData<GrnGeneType, EdgeWeight, GrnMutParams> {
    fn from(value: Grn<I, O>) -> Self {
        let nodes: Vec<Node<GrnGeneType>> = value.graph()
            .node_references()
            .map(|(i, node)| Node {
                id: i.index(),
                data: node.clone(),
            })
            .collect();

        let links: Vec<Link<EdgeWeight>> = value.graph()
            .edge_references()
            .map(|e| Link {
                source: e.source().index(),
                target: e.target().index(),
                data: EdgeWeight { weight: *e.weight() },
            })
            .collect();

        NodeLinkData {
            directed: true,
            multigraph: false,
            graph: GrnMutParams {
                mut_rate: value.mut_rate,
                mut_std: value.mut_distr.std_dev()
            },
            nodes,
            links,
        }
    }
}

impl TryFrom<NodeLinkData<GrnGeneType, EdgeWeight, GrnMutParams>> for DiGraph<GrnGeneType, f32> {
    type Error = anyhow::Error;

    fn try_from(node_link: NodeLinkData<GrnGeneType, EdgeWeight, GrnMutParams>) -> Result<Self, Self::Error> {
        if !node_link.directed {
            bail!("graph must be directed");
        }
        let mut graph = DiGraph::new();
        for node in node_link.nodes {
            graph.add_node(node.data);
        }
        for edge in node_link.links {
            graph.add_edge(
                NodeIndex::new(edge.source),
                NodeIndex::new(edge.target),
                edge.data.weight
            );
        }
        Ok(graph)
    }
}

impl<const I: usize, const O: usize> TryFrom<NodeLinkData<GrnGeneType, EdgeWeight, GrnMutParams>> for Grn<I, O> {
    type Error = anyhow::Error;

    fn try_from(value: NodeLinkData<GrnGeneType, EdgeWeight, GrnMutParams>) -> Result<Self, Self::Error> {
        if value.multigraph {
            bail!("graph can not be a multigraph");
        }
        let mut_rate = value.graph.mut_rate;
        let mut_std = value.graph.mut_std;
        Ok(Self::from_graph(
            DiGraph::try_from(value)?,
            mut_rate,
            mut_std,
        )?)
    }
}

#[derive(Serialize, Deserialize)]
pub struct GrnMutParams {
    pub(crate) mut_rate: f32,
    pub(crate) mut_std: f32
}