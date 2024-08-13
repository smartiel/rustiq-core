//! This module contains all the Pauli network synthesis algorithms.
//!
//!
//! Pauli networks are generalized parity networks. They are Clifford circuits that can be used
//! to implement a target sequence of Pauli rotations. A Clifford circuit is a Pauli network for a given sequence of
//! rotations if it can be used to implement the sequence by just inserting single qubit Pauli rotations inside it.
//!
//! Our synthesis interface can handle different settings:
//! * optimization for entangling depth or count
//! * it can preserve the rotation order or relax it (this is valid for some applications)
//! * it can restore the Clifford frame to the identity at the end of the network
//!
//! # Synthesis example
//!
//! We will synthesize a Pauli network for the sequence `["XX", "ZZ", "YY"]`.
//! ```
//! use rustiq_core::structures::{PauliSet, Metric};
//! use rustiq_core::synthesis::pauli_network::greedy_pauli_network;
//! let mut paulis = PauliSet::from_slice(&["XX".to_string(), "ZZ".to_string(), "YY".to_string()]);
//! // targeting count, not reordering, 1 attempt, no reset of the Clifford frame
//! let circuit = greedy_pauli_network(&mut paulis, &Metric::COUNT, false, 1, false, false);
//!
//! // targeting depth, not reordering, 1 attempt, no reset of the Clifford frame
//! let circuit = greedy_pauli_network(&mut paulis, &Metric::DEPTH, false, 1, false, false);
//!
//! // Same thing, relaxing rotation order
//! let circuit = greedy_pauli_network(&mut paulis, &Metric::DEPTH, true, 1, false, false);
//!
//! // Same thing, fixing the Clifford frame
//! let circuit = greedy_pauli_network(&mut paulis, &Metric::DEPTH, true, 1, false, true);
//!
//! // Same thing, performing 10 attempts and returning the best
//! let circuit = greedy_pauli_network(&mut paulis, &Metric::DEPTH, true, 10, false, true);
//! ```
pub mod chunks;
pub mod greedy_order_preserving;
pub mod greedy_pauli_network;
pub mod synthesis;

pub use synthesis::{check_circuit, greedy_pauli_network};
