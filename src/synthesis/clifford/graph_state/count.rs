use crate::routines::decoding::information_set_decoding;
/// This module contains the necessary methods to synthesize graph states
use crate::structures::clifford_circuit::{CliffordCircuit, CliffordGate};
use crate::structures::graph_state::GraphState;
use crate::structures::pauli_like::PauliLike;

#[derive(Clone, Copy)]
enum Type {
    Cz,
    Cnot,
    SCnotS,
}

type Parities = (Vec<Vec<bool>>, Vec<(Type, usize, usize)>);

fn gather_parities(circuit: &CliffordCircuit, n: usize) -> Parities {
    let mut graph_state = GraphState::new(n);
    let mut parities = Vec::new();
    let mut moves = Vec::new();
    let mut ei = vec![false; n];
    for i in 0..n {
        ei[i] = true;
        parities.push(ei.clone());
        ei[i] = false;
        moves.push((Type::Cz, 0, i));
    }
    for (index, gate) in circuit.gates.iter().enumerate() {
        graph_state.conjugate_with_gate(gate);
        match gate {
            CliffordGate::CNOT(i, j) => {
                for parity in parities.iter_mut() {
                    parity[*i] ^= parity[*j];
                }
                ei[*j] = true;
                parities.push(ei.clone());
                ei[*j] = false;
                moves.push((Type::Cz, index + 1, *j));

                parities.push(graph_state.adj[*i].clone());
                moves.push((Type::Cnot, index + 1, *i));

                graph_state.s(*i);
                parities.push(graph_state.adj[*i].clone());
                graph_state.s(*i);
                moves.push((Type::SCnotS, index + 1, *i));
            }

            CliffordGate::CZ(i, j) => {
                parities.push(graph_state.adj[*i].clone());
                moves.push((Type::Cnot, index + 1, *i));

                graph_state.s(*i);
                parities.push(graph_state.adj[*i].clone());
                graph_state.s(*i);
                moves.push((Type::SCnotS, index + 1, *i));

                parities.push(graph_state.adj[*j].clone());
                moves.push((Type::Cnot, index + 1, *j));

                graph_state.s(*j);
                parities.push(graph_state.adj[*j].clone());
                graph_state.s(*j);
                moves.push((Type::SCnotS, index + 1, *j));
            }
            _ => {}
        }
    }

    (parities, moves)
}

pub fn synthesize_graph_state_count(graph: &GraphState, niter: usize) -> CliffordCircuit {
    let mut circuit = CliffordCircuit::new(graph.n);
    for i in 0..graph.n {
        if i > 0 {
            let (parities, moves) = gather_parities(&circuit, i);
            let mut target = vec![false; i];
            target[..i].copy_from_slice(&graph.adj[i][..i]);
            let solution = information_set_decoding(&parities, &target, niter, true);
            let solution = solution.expect("Something went wrong during syndrome decoding :/");
            let mut new_circuit = CliffordCircuit::new(graph.n);
            let moves: Vec<(Type, usize, usize)> = solution
                .iter()
                .enumerate()
                .filter(|(_a, b)| **b)
                .map(|(a, _)| moves[a])
                .collect();
            for (ty, index, qbit) in moves.iter() {
                if *index == 0 {
                    match ty {
                        Type::Cnot => {
                            new_circuit.gates.push(CliffordGate::CNOT(i, *qbit));
                        }
                        Type::Cz => {
                            new_circuit.gates.push(CliffordGate::CZ(i, *qbit));
                        }
                        Type::SCnotS => {
                            new_circuit.gates.push(CliffordGate::S(*qbit));
                            new_circuit.gates.push(CliffordGate::CNOT(i, *qbit));
                            new_circuit.gates.push(CliffordGate::S(*qbit));
                        }
                    }
                }
            }
            for k in 0..circuit.gates.len() {
                new_circuit.gates.push(circuit.gates[k]);
                for (ty, index, qbit) in moves.iter() {
                    if *index == k + 1 {
                        match ty {
                            Type::Cnot => {
                                new_circuit.gates.push(CliffordGate::CNOT(i, *qbit));
                            }
                            Type::Cz => {
                                new_circuit.gates.push(CliffordGate::CZ(i, *qbit));
                            }
                            Type::SCnotS => {
                                new_circuit.gates.push(CliffordGate::S(*qbit));
                                new_circuit.gates.push(CliffordGate::CNOT(i, *qbit));
                                new_circuit.gates.push(CliffordGate::S(*qbit));
                            }
                        }
                    }
                }
            }
            circuit = new_circuit;
        }
    }

    let mut simulated = GraphState::new(graph.n);
    simulated.conjugate_with_circuit(&circuit);
    for i in 0..graph.n {
        if simulated.adj[i][i] != graph.adj[i][i] {
            circuit.gates.push(CliffordGate::S(i));
        }
    }
    circuit
}
