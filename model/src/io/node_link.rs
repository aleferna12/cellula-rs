use crate::genetics::grn::{EdgeWeight, Grn, GrnGeneType};
use rustworkx_core::petgraph::prelude::EdgeRef;
use rustworkx_core::petgraph::visit::IntoNodeReferences;
use serde::{Deserialize, Serialize};
use thiserror::Error;

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

pub trait ToNodeLink<N, E, G> {
    fn to_node_link(&self) -> NodeLinkData<N, E, G>;
}

pub trait FromNodeLink<N, E, G>: Sized {
    fn from_node_link(node_link: NodeLinkData<N, E, G>) -> Result<Self, NodeLinkError>;
}

impl<const I: usize, const O: usize> ToNodeLink<GrnGeneType, EdgeWeight, GrnMutParams> for Grn<I, O> {
    fn to_node_link(&self) -> NodeLinkData<GrnGeneType, EdgeWeight, GrnMutParams> {
        let nodes: Vec<Node<GrnGeneType>> = self.graph()
            .node_references()
            .map(|(i, node)| Node {
                id: i.index(),
                data: node.clone(),
            })
            .collect();

        let links: Vec<Link<EdgeWeight>> = self.graph()
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
                rate: self.mut_rate,
                std: self.mut_distr.std_dev()
            },
            nodes,
            links,
        }
    }
}

impl<const I: usize, const O: usize> FromNodeLink<GrnGeneType, EdgeWeight, GrnMutParams> for Grn<I, O> {
    fn from_node_link(node_link: NodeLinkData<GrnGeneType, EdgeWeight, GrnMutParams>) -> Result<Self, NodeLinkError> {
        if !node_link.multigraph {
            return Err(NodeLinkError::NotMultigraph)
        }
        if !node_link.directed {
            return Err(NodeLinkError::NotDirected)
        }
        Ok(Self::new(

        ))
    }
}

#[derive(Serialize, Deserialize)]
struct GrnMutParams {
    rate: f32,
    std: f32
}

#[derive(Error, Debug)]
pub enum NodeLinkError {
    #[error("graph must be directed")]
    NotDirected,
    #[error("graph must not be a multigraph")]
    NotMultigraph,
}