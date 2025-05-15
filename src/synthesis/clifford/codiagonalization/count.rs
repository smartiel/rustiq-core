use super::common::{make_full_rank, permute_circuit};
use crate::routines::decoding::information_set_decoding;
use crate::routines::f2_linalg::rowop;
use crate::structures::{CliffordCircuit, CliffordGate, GraphState, Metric, PauliSet};
use crate::synthesis::clifford::graph_state::synthesize_graph_state;

fn gather_parities(
    input_table: &[Vec<bool>],
    circuit: &CliffordCircuit,
    k: usize,
) -> (Vec<Vec<bool>>, Vec<(usize, usize)>) {
    let mut x_table = input_table.to_owned();
    let mut moves = Vec::new();
    let mut parities = Vec::new();
    for (i, row) in x_table.iter().enumerate().take(k) {
        parities.push(row.clone());
        moves.push((0, i));
    }
    for (index, gate) in circuit.gates.iter().enumerate() {
        if let CliffordGate::CNOT(i, j) = gate {
            rowop(&mut x_table, *i, *j);
            parities.push(x_table[*j].clone());
            moves.push((index + 1, *j));
        } else {
            panic!("This shouldn't happen");
        }
    }
    (parities, moves)
}

fn reduce_x_part(pauli_set: &PauliSet, niter: usize) -> (CliffordCircuit, Vec<usize>, GraphState) {
    let (mut circuit, row_perm, rank, (mut z_table, mut x_table)) = make_full_rank(pauli_set);
    let mut cnot_circuit = CliffordCircuit::new(pauli_set.n);
    for i in rank..x_table.len() {
        // if all(z_table[i].iter(), |x| !*x) {
        //     continue;
        // }
        let target = &x_table[i];
        let (parities, moves) = gather_parities(&x_table, &cnot_circuit, i);
        let solution = information_set_decoding(&parities, target, niter, true);
        let solution = solution.expect("Something went wrong during syndrome decoding :/");
        let moves: Vec<_> = solution
            .iter()
            .enumerate()
            .filter(|(_, b)| **b)
            .map(|(a, _)| moves[a])
            .collect();
        let mut new_circuit = CliffordCircuit::new(pauli_set.n);
        for (gindex, qbit) in moves.iter() {
            if *gindex == 0 {
                new_circuit.gates.push(CliffordGate::CNOT(*qbit, i));
            } else {
                break;
            }
        }
        for (j, gate) in cnot_circuit.gates.iter().enumerate() {
            new_circuit.gates.push(*gate);
            for (gindex, qbit) in moves.iter() {
                if *gindex == j + 1 {
                    new_circuit.gates.push(CliffordGate::CNOT(*qbit, i));
                }
            }
        }
        cnot_circuit = new_circuit;
    }
    for gate in cnot_circuit.gates.iter() {
        if let CliffordGate::CNOT(i, j) = gate {
            rowop(&mut z_table, *j, *i);
            rowop(&mut x_table, *i, *j);
        }
    }
    let permuted_circuit = permute_circuit(&cnot_circuit, &row_perm);

    let mut graph_state = GraphState::new(rank);
    for col in 0..rank {
        for (row_index, row) in z_table.iter().enumerate().take(rank) {
            graph_state.adj[row_index][col] = row[col];
        }
    }
    circuit.extend_with(&permuted_circuit);
    (circuit, row_perm, graph_state)
}

pub fn codiagonalize_count(pauli_set: &PauliSet, niter: usize) -> CliffordCircuit {
    let (mut circuit, perm, graph) = reduce_x_part(pauli_set, niter);
    let gs_synth = synthesize_graph_state(&graph, &Metric::COUNT, niter);
    let gs_synth = permute_circuit(&gs_synth, &perm);
    circuit.extend_with(&gs_synth.dagger());
    for bit in perm.iter().take(graph.n) {
        circuit.gates.push(CliffordGate::H(*bit));
    }
    circuit
}

#[cfg(test)]
mod codiag_count_tests {

    use super::*;
    use crate::structures::pauli_like::PauliLike;
    #[test]
    fn test_synth() {
        let mut pset = PauliSet::new(4);
        pset.insert("ZZII", false);
        pset.insert("IIXI", false);
        pset.insert("IIIZ", false);
        codiagonalize_count(&pset, 1);
    }

    #[test]
    fn test_cnot_part() {
        let mut pset = PauliSet::new(7);
        pset.insert("XIIXIXX", false);
        pset.insert("IXIXXII", false);
        pset.insert("IIXXIII", false);
        let circuit = codiagonalize_count(&pset, 1);
        println!("{:?}", circuit);
    }
    use rand::Rng;
    fn random_instance(n: usize, m: usize) -> PauliSet {
        let mut rng = rand::rng();
        let mut pset = PauliSet::new(n);
        for _ in 0..m {
            let mut vec: Vec<bool> = vec![false; 2 * n];
            for b in vec.iter_mut().take(n) {
                *b = rng.random::<bool>();
            }

            pset.insert_vec_bool(&vec, false);
        }
        for _ in 0..n * n {
            let i = rng.random::<u64>() as usize % n;
            loop {
                let j = rng.random::<u64>() as usize % n;
                if j != i {
                    pset.cnot(i, j);
                    let g2 = rng.random::<bool>();
                    if g2 {
                        pset.h(j);
                    } else {
                        pset.s(j);
                    }
                    break;
                }
            }
            let g1 = rng.random::<bool>();
            if g1 {
                pset.h(i);
            } else {
                pset.s(i);
            }
        }
        pset
    }

    #[test]
    fn test_thin() {
        for _ in 0..10 {
            let instance = random_instance(50, 20);
            let mut copy_instance = instance.clone();
            let circuit = codiagonalize_count(&instance, 1);
            copy_instance.conjugate_with_circuit(&circuit);
            for i in 0..instance.len() {
                let (_, vec) = copy_instance.get_as_vec_bool(i);
                assert!(vec[..instance.n].iter().all(|b| !*b));
            }
        }
    }
    #[test]
    fn test_thick() {
        for _ in 0..10 {
            let instance = random_instance(20, 50);
            let mut copy_instance = instance.clone();
            let circuit = codiagonalize_count(&instance, 1);
            copy_instance.conjugate_with_circuit(&circuit);
            for i in 0..instance.len() {
                let (_, vec) = copy_instance.get_as_vec_bool(i);
                assert!(vec[..instance.n].iter().all(|b| !*b));
            }
        }
    }
}
