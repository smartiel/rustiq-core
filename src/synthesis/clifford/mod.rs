//! This module contains all the Clifford circuit synthesis algorithms.
//!
//! * [Clifford and Clifford isometry synthesis ](isometry::isometry_synthesis)
//! * [Pauli operators codiagonalization](codiagonalization::codiagonalize)
//! * [Graph states](graph_state::synthesize_graph_state) and [stabilizer states synthesis](graph_state::synthesize_stabilizer_state)
//!
//! # Clifford synthesis examples
//!
//! `Tableau` synthesis is handled as follows:
//!
//! ```
//! use rustiq_core::structures::{Tableau, Metric};
//! use rustiq_core::synthesis::clifford::isometry::isometry_synthesis;
//! // Allocating a random Tableau over 5 qubits.
//! let my_tab = Tableau::random(5);
//! let circuit = isometry_synthesis(&my_tab.to_isometry(), &Metric::COUNT, 10);
//!
//! println!("{:?}", circuit.entangling_depth());
//! println!("{:?}", circuit.entangling_count());
//! ```
//!
//! The synthesis function is also able to target entangling depth instead of count:
//!```
//! use rustiq_core::structures::{Tableau, Metric};
//! use rustiq_core::synthesis::clifford::isometry::isometry_synthesis;
//! // Allocating a random Tableau over 5 qubits.
//! let my_tab = Tableau::random(5);
//! let circuit = isometry_synthesis(&my_tab.to_isometry(), &Metric::DEPTH, 0);
//!
//! println!("{:?}", circuit.entangling_depth());
//! println!("{:?}", circuit.entangling_count());
//! ```
//!
//! The same synthesis function can also be used to optimize a Clifford circuit:
//!
//! ```
//! use rustiq_core::structures::{Tableau, Metric, CliffordCircuit, CliffordGate};
//! use rustiq_core::synthesis::clifford::isometry::isometry_synthesis;
//! // Creating a sub-optimal circuit:
//! let mut my_circuit = CliffordCircuit::new(2);
//! my_circuit.gates.push(CliffordGate::CNOT(0, 1));
//! my_circuit.gates.push(CliffordGate::CNOT(0, 1));
//!
//! // Generating the corresponding `Tableau`
//! let my_tab = Tableau::from_circuit(&my_circuit);
//!
//! // Calling the synthesis algorithm
//! let circuit = isometry_synthesis(&my_tab.to_isometry(), &Metric::COUNT, 1);
//!
//! assert_eq!(circuit.cnot_count(), 0);
//! ```
//!
//! # Pauli codiagonalization example
//!
//! This algorithms acts on a set of pairwise commuting Pauli operators.
//! In this example, we will codiagonalize operators `XXX, ZXY, YXZ`.
//!
//! Similarly to all the other synthesis algorithms, we can either target entangling count or depth.
//! ```
//! use rustiq_core::structures::{PauliSet, Metric, PauliLike};
//! use rustiq_core::synthesis::clifford::codiagonalization::codiagonalize;
//!
//! // Storing our Pauli operators in a `PauliSet` struct
//! let mut paulis = PauliSet::from_slice(&["XXX".to_string(), "ZXY".to_string(), "YXZ".to_string()]);
//!
//! // Depth call
//! let circuit = codiagonalize(&mut paulis.clone(), &Metric::DEPTH, 0);
//! // Count call: the last parameter can be increase to improve output quality (slowing down the algorithm)
//! let circuit = codiagonalize(&mut paulis.clone(), &Metric::COUNT, 10);
//!
//! paulis.conjugate_with_circuit(&circuit);
//!
//! println!("{}", paulis);
//!
//! ```
//!

pub mod codiagonalization;
pub mod graph_state;
pub mod isometry;
