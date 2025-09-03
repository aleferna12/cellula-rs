use rustworkx_core::petgraph::prelude::EdgeRef;
use rustworkx_core::petgraph::visit::IntoNodeReferences;
use serde::Serialize;
use crate::genetics::grn::{EdgeWeight, Grn, GrnGeneType};

#[derive(Serialize)]
pub struct Node<N> {
    pub id: usize,
    #[serde(flatten)]
    pub data: N,
}

#[derive(Serialize)]
pub struct Link<E> {
    pub source: usize,
    pub target: usize,
    #[serde(flatten)]
    pub data: E,
}

#[derive(Serialize)]
pub struct NodeLinkData<N, E> {
    pub directed: bool,
    pub multigraph: bool,
    pub graph: serde_json::Map<String, serde_json::Value>,
    pub nodes: Vec<Node<N>>,
    pub links: Vec<Link<E>>,
}

pub trait ToNodeLink {
    type Node;
    type Edge;
    fn to_node_link(&self) -> NodeLinkData<Self::Node, Self::Edge>;
}

impl<const I: usize, const O: usize> ToNodeLink for Grn<I, O> {
    type Node = GrnGeneType;
    type Edge = EdgeWeight;

    fn to_node_link(&self) -> NodeLinkData<Self::Node, Self::Edge> {
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
            graph: serde_json::Map::new(),
            nodes,
            links,
        }
    }
}