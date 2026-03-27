use pyo3::prelude::*;
use std::collections::HashMap;

#[pymodule]
pub mod neighbor {
    #[pymodule_export]
    use super::neighbor_set_index;
}

/// Returns a Vec<usize> of the same length as nodes, where each element
/// is the graph ID that the corresponding node belongs to.
///
/// Graph IDs are assigned in order of first appearance (0, 1, 2, ...).
#[pyfunction]
pub fn neighbor_set_index(spins: Vec<String>, neighbors: Vec<Vec<String>>) -> Vec<usize> {
    find_subgraphs(&spins, &neighbors)
}

fn find_subgraphs<T>(nodes: &[T], edges: &[Vec<T>]) -> Vec<usize>
where
    T: Eq + std::hash::Hash + Clone,
{
    let node_to_idx: HashMap<T, usize> = nodes
        .iter()
        .cloned()
        .enumerate()
        .map(|(i, n)| (n, i))
        .collect();

    let n = nodes.len();
    let mut parent: Vec<usize> = (0..n).collect();
    let mut rank: Vec<usize> = vec![0; n];

    for (i, neighbors) in edges.iter().enumerate() {
        for neighbor in neighbors {
            let j = node_to_idx[neighbor];
            union(&mut parent, &mut rank, i, j);
        }
    }

    // Assign a stable graph ID to each unique root, in order of first appearance
    let mut root_to_graph_id: HashMap<usize, usize> = HashMap::new();
    let mut next_id = 0;

    nodes
        .iter()
        .enumerate()
        .map(|(i, _)| {
            let root = find(&mut parent, i);
            *root_to_graph_id.entry(root).or_insert_with(|| {
                let id = next_id;
                next_id += 1;
                id
            })
        })
        .collect()
}

fn find(parent: &mut [usize], x: usize) -> usize {
    let mut root = x;
    while parent[root] != root {
        root = parent[root];
    }
    // Path compression
    let mut curr = x;
    while parent[curr] != root {
        let next = parent[curr];
        parent[curr] = root;
        curr = next;
    }
    root
}

fn union(parent: &mut [usize], rank: &mut [usize], a: usize, b: usize) {
    let ra = find(parent, a);
    let rb = find(parent, b);
    if ra == rb {
        return;
    }
    let (ra, rb) = if rank[ra] < rank[rb] { (rb, ra) } else { (ra, rb) };
    parent[rb] = ra;
    if rank[ra] == rank[rb] {
        rank[ra] += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_two_components() {
        let nodes = vec![10, 20, 30, 40, 50];
        let edges = vec![vec![20], vec![10, 30], vec![20], vec![50], vec![40]];
        assert_eq!(find_subgraphs(&nodes, &edges), vec![0, 0, 0, 1, 1]);
    }

    #[test]
    fn test_all_connected() {
        let nodes = vec![1, 2, 3, 4];
        let edges = vec![vec![2], vec![1, 3], vec![2, 4], vec![3]];
        assert_eq!(find_subgraphs(&nodes, &edges), vec![0, 0, 0, 0]);
    }

    #[test]
    fn test_no_edges() {
        let nodes = vec![7, 8, 9];
        let edges = vec![vec![], vec![], vec![]];
        assert_eq!(find_subgraphs(&nodes, &edges), vec![0, 1, 2]);
    }

    #[test]
    fn test_single_node() {
        let nodes = vec![42];
        let edges = vec![vec![]];
        assert_eq!(find_subgraphs(&nodes, &edges), vec![0]);
    }

    #[test]
    fn test_three_components() {
        let nodes = vec![5, 10, 15, 20, 25, 30];
        let edges = vec![vec![10], vec![5], vec![20], vec![15], vec![], vec![]];
        assert_eq!(find_subgraphs(&nodes, &edges), vec![0, 0, 1, 1, 2, 3]);
    }
}