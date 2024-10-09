#[allow(unused_imports)]
use crate::structures::{CliffordCircuit, Metric, PauliDag, PauliLike, PauliSet};
/// This module implements a data structure encoding the commutation DAG of a sequence of Pauli operators.

pub fn pauli_network_synthesis_no_permutation(
    axes: &mut PauliSet,
    metric: &Metric,
    skip_sort: bool,
) -> CliffordCircuit {
    let mut dag = PauliDag::from_pauli_set(axes.clone());
    let mut circuit = CliffordCircuit::new(dag.pauli_set.n);
    while !dag.front_layer.is_empty() {
        circuit.extend_with(&dag.single_step_synthesis(metric, skip_sort));
    }
    circuit
}

#[cfg(test)]
mod greedy_synthesis_tests {
    use super::*;
    use crate::structures::clifford_circuit::CliffordGate;
    use std::collections::HashSet;
    fn check_circuit(input: &[String], circuit: &CliffordCircuit) -> bool {
        let mut hit_map: HashSet<usize> = HashSet::new();
        let mut bucket = PauliSet::from_slice(input);
        for i in 0..bucket.len() {
            if bucket.support_size(i) == 1 {
                hit_map.insert(i);
            }
        }
        for gate in circuit.gates.iter() {
            match gate {
                CliffordGate::SqrtX(i) => {
                    bucket.sqrt_x(*i);
                }
                CliffordGate::SqrtXd(i) => {
                    bucket.sqrt_xd(*i);
                }
                CliffordGate::S(i) => {
                    bucket.s(*i);
                }
                CliffordGate::Sd(i) => {
                    bucket.sd(*i);
                }
                CliffordGate::H(i) => {
                    bucket.h(*i);
                }
                CliffordGate::CNOT(i, j) => {
                    bucket.cnot(*i, *j);
                }
                _ => (),
            }
            for i in 0..bucket.len() {
                if bucket.support_size(i) == 1 {
                    hit_map.insert(i);
                }
            }
        }
        println!("Synthesized {} operators", hit_map.len());
        println!("{:?}", bucket);
        return hit_map.len() == input.len();
    }
    #[test]
    fn count_synthesis() {
        let axes = vec!["XX".to_owned(), "ZZ".to_owned(), "YY".to_owned()];
        let mut input = PauliSet::from_slice(&axes);
        let circuit = pauli_network_synthesis_no_permutation(&mut input, &Metric::COUNT, false);
        println!("{circuit:?}");
        assert!(check_circuit(&axes, &circuit))
    }
    #[test]
    fn depth_synthesis() {
        let axes = vec!["XX".to_owned(), "ZZ".to_owned(), "YY".to_owned()];
        let mut input = PauliSet::from_slice(&axes);
        let circuit = pauli_network_synthesis_no_permutation(&mut input, &Metric::DEPTH, false);
        println!("{circuit:?}");
        assert!(check_circuit(&axes, &circuit))
    }

    #[test]
    fn count_synthesis_ss() {
        let axes = vec!["XX".to_owned(), "ZZ".to_owned(), "YY".to_owned()];
        let mut input = PauliSet::from_slice(&axes);
        let circuit = pauli_network_synthesis_no_permutation(&mut input, &Metric::COUNT, true);
        println!("{circuit:?}");
        assert!(check_circuit(&axes, &circuit))
    }
    #[test]
    fn depth_synthesis_ss() {
        let axes = vec!["XX".to_owned(), "ZZ".to_owned(), "YY".to_owned()];
        let mut input = PauliSet::from_slice(&axes);
        let circuit = pauli_network_synthesis_no_permutation(&mut input, &Metric::DEPTH, true);
        println!("{circuit:?}");
        assert!(check_circuit(&axes, &circuit))
    }
}
