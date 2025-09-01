use serde::Serialize;

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
    pub graph: serde_json::Value, // empty object {}
    pub nodes: Vec<Node<N>>,
    pub links: Vec<Link<E>>,
}

pub trait ToNodeLink {
    type Node;
    type Edge;
    fn to_node_link(&self) -> NodeLinkData<Self::Node, Self::Edge>;
}