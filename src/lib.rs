//! `rustiq-core` is a quantum circuit synthesis library.
//!
//! It contains data structures representing various mathematical objects together with synthesis methods able to generate
//! efficient circuits targeting 2-qubit gate count or depth.
//!
//! # Data structures
//!
//! * [CliffordCircuit](crate::structures::CliffordCircuit) and [CliffordGate](crate::structures::CliffordGate) - Data structures used to represent Clifford gates and circuits.
//! * [PauliSet](crate::structures::PauliSet) - A data structure used to store a list of Pauli operators. The list can be efficiently conjugated by `CliffordGate` or `CliffordCircuit` objects.
//! * [Tableau](crate::structures::Tableau) - A Clifford Tableau implementation built on top of `PauliSet`.
//! * [IsometryTableau](crate::structures::IsometryTableau) - An extension of the Tableau data structure that describes a Clifford operator applied to a partially stabilized input.
//! * [GraphState](crate::structures::GraphState) - A data structure representing a graph state
//!
//! # Synthesis algorithms
//!
//! All synthesis algorithms are located in the [synthesis] submodule.
//!
//! * Synthesis of Clifford operators are handled by the [isometry_synthesis](crate::synthesis::clifford::isometry::isometry_synthesis) method.
//!   This method can handle synthesis of `IsometryTableau` objects into `CliffordCircuit`.
//! * Graph states and stabilizer states can be synthesized using [synthesize_graph_state](crate::synthesis::clifford::graph_state::synthesize_graph_state)
//!   and [synthesize_stabilizer_state](crate::synthesis::clifford::graph_state::synthesize_stabilizer_state).
//! * Method [codiagonalize](crate::synthesis::clifford::codiagonalization::codiagonalize) can be used to produce a
//!   `CliffordCircuit` that codiagonalizes a given set of pairwise commuting Pauli operators.
//! * Method [greedy_pauli_network](crate::synthesis::pauli_network::greedy_pauli_network()) can be used to implement a sequence of Pauli rotations.
//!

pub mod routines;
pub mod structures;
pub mod synthesis;
