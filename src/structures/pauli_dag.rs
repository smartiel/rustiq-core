use super::{CliffordCircuit, Metric, PauliLike, PauliSet};
use crate::synthesis::pauli_network::greedy_pauli_network::single_synthesis_step;
use petgraph::prelude::*;

pub type Dag = DiGraph<usize, ()>;

/// Constructs an anti-commutation Dag from a set of operators
pub fn build_dag_from_pauli_set(pauli_set: &PauliSet) -> Dag {
    let mut dag = Dag::new();
    let node_indices: Vec<NodeIndex> = (0..pauli_set.len()).map(|i| dag.add_node(i)).collect();

    for i in 0..pauli_set.len() {
        for j in 0..i {
            if !pauli_set.commute(i, j) {
                dag.add_edge(node_indices[i], node_indices[j], ());
            }
        }
    }
    dag
}

/// Computes the list of operators that can be synthesized
pub fn get_front_layer(dag: &Dag) -> Vec<NodeIndex> {
    return dag
        .node_indices()
        .filter(|node| dag.neighbors(*node).collect::<Vec<_>>().is_empty())
        .collect();
}

pub struct PauliDag {
    /// A global set containing all the operators
    pub pauli_set: PauliSet,
    /// The dag structure
    pub dag: Dag,
    /// The front layer of (unprocessed) DAG nodes
    /// (corresponding Pauli operators pairwise commute)
    pub front_nodes: Vec<NodeIndex>,
    /// Stores the number of (unprocessed) predecessors for each node,
    /// indexed by node id
    pub in_degree: Vec<usize>,
}

impl PauliDag {
    /// Constructs a PauliDag from a PauliSet
    pub fn from_pauli_set(pauli_set: PauliSet) -> Self {
        let dag = build_dag_from_pauli_set(&pauli_set);
        let num_nodes = dag.node_count();

        let mut in_degree: Vec<usize> = vec![0; num_nodes];
        let mut front_nodes: Vec<NodeIndex> = Vec::new();
        for node_index in dag.node_indices() {
            let node_in_degree = dag.neighbors_directed(node_index, Incoming).count();
            in_degree[node_index.index()] = node_in_degree;
            if node_in_degree == 0 {
                front_nodes.push(node_index);
            }
        }

        Self {
            pauli_set,
            dag,
            front_nodes,
            in_degree,
        }
    }

    /// Constructs a PauliDag from a slice of axes
    pub fn from_slice(axes: &[String]) -> Self {
        Self::from_pauli_set(PauliSet::from_slice(axes))
    }

    /// Checks if the Pauli corresponding to `node_index` is synthesized,
    /// that is, if its support is of size <= 1
    fn is_synthesized(&self, node_index: NodeIndex) -> bool {
        self.pauli_set
            .support_size(*self.dag.node_weight(node_index).unwrap())
            <= 1
    }

    /// Remove synthesized operations from the front layer and update
    /// the front layer accordingly.
    pub(crate) fn update_front_nodes(&mut self) {
        let mut unprocessed = self.front_nodes.clone();
        self.front_nodes = Vec::new();
        // For some reason this is less performant than cloning:
        // let mut unprocessed: Vec<NodeIndex> = Vec::new();
        // std::mem::swap(&mut self.front_nodes, &mut unprocessed);

        while let Some(node_index) = unprocessed.pop() {
            if !self.is_synthesized(node_index) {
                self.front_nodes.push(node_index);
            } else {
                // the node can be removed, check which of its successors are now
                // front nodes
                for successor in self.dag.neighbors_directed(node_index, Outgoing) {
                    self.in_degree[successor.index()] -= 1;
                    if self.in_degree[successor.index()] == 0 {
                        unprocessed.push(successor);
                    }
                }
            }
        }
    }

    /// Returns true if fully processed
    pub fn fully_processed(&self) -> bool {
        self.front_nodes.is_empty()
    }

    /// Performs a single synthesis step
    pub fn single_step_synthesis(
        &mut self,
        metric: &Metric,
        skip_sort: bool,
        synthesized_circuit: &mut CliffordCircuit,
    ) {
        // Creating a fresh PauliSet from the nodes in the front layer
        let mut front_layer: PauliSet = PauliSet::new(self.pauli_set.n);
        for index in &self.front_nodes {
            let (phase, pstring) = self.pauli_set.get(*self.dag.node_weight(*index).unwrap());
            front_layer.insert(&pstring, phase);
        }

        if !skip_sort {
            front_layer.support_size_sort();
        }

        let circuit = single_synthesis_step(&mut front_layer, metric);

        // Updating the global set of operators
        self.pauli_set.conjugate_with_circuit(&circuit);

        synthesized_circuit.extend_with(&circuit);

        // Updating the front layer
        self.update_front_nodes();
    }
}
